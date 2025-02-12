#!/usr/bin/env python3
import os
from glob import glob
import argparse
from collections import defaultdict
import warnings
from mskpy.image.core import rebin

import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import log as astropy_log
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS, FITSFixedWarning
import sep
from sbpy.data import Ephem, Orbit, TimeScaleWarning, QueryError

from calviacat import RefCat2
from mskpy import gcentroid
from quick_look import (
    rho_labels,
    rho_km,
    rho_arcsec,
    skip_targets,
    rename_target,
    locations,
    color_corrections,
    assumed_gmr,
    filters,
    setup_logger,
)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("sources", nargs="*", help="measure these files (or targets with --target)")
parser.add_argument("--target", dest="mode", default="file", action="store_const", const="target", help="process targets, rather than files",)
parser.add_argument(
    "-f", "--force", action="store_true", help="re-measure and overwrite existing data"
)
parser.add_argument(
    "-v", dest="level", action="store_const", const="DEBUG", help="verbose mode"
)
parser.add_argument(
    "-q", dest="level", action="store_const", const="WARNING", help="quiet mode"
)
args = parser.parse_args()

if args.level != "DEBUG":
    warnings.simplefilter("ignore", FITSFixedWarning)
    warnings.simplefilter("ignore", RuntimeWarning)

logger = setup_logger(__file__, args.level)
logger.info("Start")

# prevents logging of fitter warnings from calviacat
astropy_log.setLevel("ERROR")

if not os.path.exists("zeropoints"):
    os.system("mkdir zeropoints")

if not os.path.exists("backgrounds"):
    os.system("mkdir backgrounds")

# skip_objects = ['C/2021 C4']  # not in horizons
files = glob("e91/*/*fz")
if args.mode == "file" and len(args.sources) > 0:
    files = args.sources
    # validate input
    for f in files:
        if not f.startswith("e91/"):
            raise ValueError(f"{f} must start with e91/")

if os.path.exists("phot-skip.list"):
    skip = ascii.read("phot-skip.list")
else:
    skip = Table(names=("file", "reason"), dtype=("U64", "U128"))


bgopts = dict(bw=256, bh=256, fw=5, fh=5)


def set_table_formats(tab):
    tab["fracday"].format = "{:.5f}"
    tab["tmtp"].format = "{:.3f}"
    tab["rh"].format = "{:.3f}"
    tab["delta"].format = "{:.3f}"
    tab["phase"].format = "{:.3f}"
    tab["pixel scale"].format = "{:.3f}"
    tab["zp"].format = "{:.4f}"
    tab["color cor"].format = "{:.4f}"
    tab["zp err"].format = "{:.4f}"
    tab["exptime"].format = "{:.1f}"
    tab["airmass"].format = "{:.3f}"
    tab["seeing"].format = "{:.2f}"
    tab["rho10k"].format = "{:.2f}"
    tab["cx"].format = "{:.1f}"
    tab["cy"].format = "{:.1f}"
    tab["dc"].format = "{:.1f}"
    tab["bgarea"].format = "{:.0f}"
    tab["bg"].format = "{:.5g}"
    tab["bgsig"].format = "{:.5g}"
    for r in rho_labels:
        tab["flux{}".format(r)].format = "{:.5g}"
        tab["m{}".format(r)].format = "{:.3f}"
        tab["merr{}".format(r)].format = "{:.3f}"
    return tab


if os.path.exists("phot.txt"):
    tab = set_table_formats(ascii.read("phot.txt"))
else:
    tab = {"file": []}

# ephemeris offset
offset = defaultdict(lambda: {"ra": 0 * u.arcsec, "dec": 0 * u.arcsec})
offset["141P"] = {"ra": 1.4 * u.arcsec, "dec": -0.9 * u.arcsec}

# centroiding box (box sizes are binned pixels)
centroid_options = defaultdict(lambda: {"binning": 1, "boxes": [17, 11]})
#centroid_options["C/2014 UN271"] = {"binning": 2, "boxes": None}
# centroid_options['7P'] = {
#     'binning': 1,
#     'boxes': [9, 5]
# }
# centroid_boxes = defaultdict(lambda: [17, 11])
# centroid_boxes['2014 UN271'] = [21]

new_objects = set()

rc2 = RefCat2("cat.db", min_matches=10, logger=logger)

for f in files:
    basename = os.path.basename(f)[:-8]
    if f in tab["file"] and basename in skip["file"]:
        # this file was added to the skip list after photometry; remove it
        if isinstance(tab, dict):
            i = tab["file"].index(f)
            del tab["file"][i]
        else:
            i = tab["file"] != f
            tab = tab[i]

    if basename in skip["file"]:
        continue

    # file mode? file already processed? force reprocessing not enabled?
    if args.mode == "file" and f in tab["file"] and not args.force:
        logger.debug("file mode, but file already processed and --force not enabled: %s", f)
        continue

    bgf = "backgrounds/{}.fits.gz".format(basename)

    try:
        hdu = fits.open(f)
    except OSError:
        logger.exception("Error opening %s.", basename)
        skip.add_row((basename, "Error opening file"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue
    h = hdu["sci"].header

    target = rename_target.get(h["OBJECT"], h["OBJECT"])
    if target in skip_targets:
        logger.debug("target in skip_targets: %s", f)
        continue

    if args.mode == "target":
        if target not in args.sources:
            logger.debug("target mode and this target was not requested: %s", f)
            continue
        elif f in tab["file"] and not args.force:
            logger.debug("target mode but this file already processed and --force not enabled: %s", f)
            continue

    # at this point if the file is in the table, we need to remove it
    if f in tab["file"]:
        tab = tab[tab["file"] != f]

    bpm = hdu["bpm"].data != 0
    im = hdu["sci"].data + 0
    # do not use gain, already in e-
    err = np.sqrt((im + h["BIASLVL"]) + h["RDNOISE"] ** 2)

    if "cat" not in hdu:
        logger.exception("%s missing photometry catalog.", basename)
        skip.add_row((basename, "Missing photometry catalog"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    phot = Table(hdu["cat"].data)
    if "ra" not in phot.colnames:
        logger.error("Adding %s to skip file list: WCS probably missing.", basename)
        skip.add_row((basename, "WCS probably missing"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    phot = phot[phot["flux"] > 0]

    new_objects.add(target)
    ps = h["PIXSCALE"] * u.arcsec

    # first determine if target is in the image
    midtime = Time(h["DATE-OBS"]) + h["EXPTIME"] / 2 * u.s
    opts = dict(epochs=midtime, cache=True, location=locations[basename[:3]])
    eph = None
    if target.startswith("C/") or target.startswith("P/") or target.endswith("P"):
        opts2 = dict(
            id_type="designation",
            closest_apparition=True,
            no_fragments=True,
        )
    else:
        opts2 = dict(id_type="smallbody")

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", (UserWarning, TimeScaleWarning))
            eph = Ephem.from_horizons(target, **opts, **opts2)
            orb = Orbit.from_horizons(
                target,
                epochs=opts["epochs"],
                cache=opts["cache"],
                center="Sol",
                **opts2,
            )
    except QueryError:
        logger.error(
            "Ephemeris request error.  Adding to skip list.  Add a name translation to lco_util and restore?"
        )
        skip.add_row((basename, "ephemeris request error, possible issue with name"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    eph["rh"] = [np.sign(row["rdot"]) * row["rh"] for row in eph]
    Tp = orb["Tp"].jd

    # target location
    wcs = WCS(h)
    c = wcs.all_world2pix(
        eph["ra"] + offset[target]["ra"], eph["dec"] + offset[target]["dec"], 0
    )
    if any(np.isnan(c).ravel()):
        logger.error(
            "Adding %s to skip file list: WCS transformation error."
            "  Ephemeris: %s, %s",
            basename,
            eph["ra"] + offset[target]["ra"],
            eph["dec"] + offset[target]["dec"],
        )
        skip.add_row((basename, "WCS transformation error"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    gx, gy = np.array(c).ravel()

    if min(gx, gy) < 0 or gx > im.shape[1] or gy > im.shape[0]:
        logger.error("Adding %s to skip file list: Target not in frame.", basename)
        skip.add_row((basename, "Target not in frame"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    # source mask
    if not os.path.exists(bgf):
        # quick bg filter
        mask = bpm.copy()
        bkg = sep.Background(im, mask=mask, **bgopts)
        _im = im - bkg.back()

        mask += nd.grey_closing(_im, 5) > bkg.globalrms * 2
    else:
        mask = fits.getdata(bgf, "mask").astype(bool)

    source_mask = mask
    del mask

    # calibrate
    lco = SkyCoord(phot["ra"], phot["dec"], unit="deg")

    if len(rc2.search(lco)[0]) < 100:
        rc2.fetch_field(lco)

    try:
        objids, distances = rc2.xmatch(lco)
    except TypeError:
        logger.error("Adding %s to skip file list: Error in xmatch.", basename)
        skip.add_row((basename, "Error in catalog xmatch"))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    m_inst = -2.5 * np.log10(phot["flux"])
    m_err = phot["fluxerr"] / phot["flux"] * 1.0857
    catfilt = filters.get(h["FILTER"])
    if catfilt is None:
        msg = "No calibration path for filter " + h["FILTER"]
        logger.exception(msg)
        skip.add_row((basename, msg))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue

    C = color_corrections.get(h["telid"], {}).get(h["filter"], None)
    try:
        zp, C, zp_unc, m, gmr, gmi = rc2.cal_color(
            objids, m_inst, catfilt, "g-r", C=C, mlim=[13, 18], gmi_lim=[0.2, 3.0]
        )
    except Exception as e:
        logger.exception("Calibration error with %s", basename)
        skip.add_row((basename, "Calibration error."))
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        continue
    n_cal = (~m.mask).sum()

    # plot
    fig = plt.figure(1)
    fig.clear()
    ax = fig.gca()

    ax.scatter(gmr, m - m_inst, marker=".", color="k")

    x = np.linspace(0, 1.5)
    label = "{:.4f} + {:.4f} ($g-r$)".format(zp, C)
    ax.plot(x, C * x + zp, "r-", label=label)

    plt.setp(ax, xlabel="$g-r$ (mag)", ylabel=r"$m-m_{\rm inst}$ (mag)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("zeropoints/{}.png".format(basename))

    # 2D background
    if not os.path.exists(bgf):
        bkg = sep.Background(im, mask=source_mask, **bgopts)
        bghdu = fits.HDUList()
        bghdu.append(fits.PrimaryHDU(bkg.back().astype(np.float32)))
        bghdu.append(fits.ImageHDU(source_mask.astype(np.uint8), name="mask"))
        bghdu[0].header["median"] = (bkg.globalback, "global median")
        bghdu[0].header["medrms"] = (bkg.globalrms, "global rms")
        bghdu.writeto(bgf, overwrite=True)

    bgim, bgh = fits.getdata(bgf, header=True)
    im -= bgim
    bg = 0
    bgsig = bgh["medrms"]
    bgarea = np.sum(~source_mask)

    binning = centroid_options[target]["binning"]
    boxes = centroid_options[target]["boxes"]
    if boxes is None:
        cy = gy
        cx = gx
        dc = np.ma.masked
    else:
        binned = rebin(im, -binning, trim=True)
        gyb = (gy + 0.5) / binning - 0.5
        gxb = (gx + 0.5) / binning - 0.5
        cyb, cxb = gyb, gxb
        for box in boxes:
            try:
                cyb, cxb = gcentroid(binned, (cyb, cxb), box=box)
                dc = np.hypot(cxb - gxb, cyb - gyb) * binning
            except:
                dc = np.ma.masked

        cy = (cyb + 0.5) * binning - 0.5
        cx = (cxb + 0.5) * binning - 0.5

    # for box in centroid_boxes.get(target):
    #     try:
    #         cy, cx = gcentroid(im, (cy, cx), box=box)
    #         dc = np.hypot(cx - gx, cy - gy)
    #     except:
    #         dc = np.ma.masked

    # background annulus
    # bgap = pu.CircularAnnulus((cx, cy), 30, 50)
    # mask = bgap.to_mask(method='center')
    # b = mask.multiply(im)
    # mms = sigma_clipped_stats(b[b != 0])
    # bg = mms[1]
    # #bg = 0
    # bgsig = mms[2]
    # bgarea = bgap.area

    rho10k = float(
        (1e4 * u.km / (725 * u.km / u.arcsec / u.au) / eph["delta"])
        .to(u.arcsec)
        .value[0]
    )
    rho = [r * u.arcsec for r in rho_arcsec]
    rho.extend([r / 725 / float(eph["delta"].value[0]) * u.arcsec for r in rho_km])
    rho = u.Quantity(rho, u.arcsec)
    rap = (rho / ps).value
    area = np.pi * rap**2

    flux, err, flags = sep.sum_circle(im, [cx], [cy], rap, err=err)
    flux -= area * bg
    err2 = err**2 + area * bgsig**2 * (1 + area / bgarea)

    mc = -2.5 * np.log10(flux) + zp + C * assumed_gmr[target]
    mc_err = np.sqrt(err2 / (flux * 1.0857) ** 2 + zp_unc**2)

    t = Time(midtime)
    fracday = t.mjd % 1
    tmtp = t.jd - Tp

    row = (
        [
            f,
            target,
            midtime.iso[:-4],
            fracday,
            tmtp[0],
            float(eph["ra"].value[0]),
            float(eph["dec"].value[0]),
            float(eph["r"].value[0]),
            float(eph["delta"].value[0]),
            float(eph["alpha"].value[0]),
            ps.value,
            h["FILTER"],
            catfilt,
            n_cal,
            zp,
            C,
            zp_unc,
            h["EXPTIME"],
            h["AIRMASS"],
            h["L1FWHM"],
            rho10k,
            cx,
            cy,
            dc,
            bgarea,
            bg,
            bgsig,
        ]
        + list(flux)
        + list(mc)
        + list(mc_err)
    )

    if isinstance(tab, dict):
        names = [
            "file",
            "target",
            "date",
            "fracday",
            "tmtp",
            "ra",
            "dec",
            "rh",
            "delta",
            "phase",
            "pixel scale",
            "filter",
            "catalog filter",
            "N cal",
            "zp",
            "color cor",
            "zp err",
            "exptime",
            "airmass",
            "seeing",
            "rho10k",
            "cx",
            "cy",
            "dc",
            "bgarea",
            "bg",
            "bgsig",
        ]
        names.extend(["flux{}".format(r) for r in rho_labels])
        names.extend(["m{}".format(r) for r in rho_labels])
        names.extend(["merr{}".format(r) for r in rho_labels])
        tab = Table(rows=[row], names=names, masked=True)
        tab.meta["comments"] = (
            [
                Time.now().iso,
                "Filters calibrated to PS1: " + str(filters),
                "Color corrected assuming g-r = ",
            ]
            + ["  {}: {}".format(k, v) for k, v in assumed_gmr.items()]
            + [
                'm2, m5, etc. use 2", 5", etc radius apertures',
                "m10k, etc. use 10,000 km etc. radius apertures",
            ]
        )
        tab = set_table_formats(tab)
    else:
        tab.add_row(row)

    tab.sort(("target", "filter", "date"))
    tab.write("phot.txt", format="ascii.fixed_width_two_line", overwrite=True)

if len(new_objects) > 0:
    logger.info("New observations of %s", new_objects)
else:
    logger.info("No new observations")

logger.info("Finish")
