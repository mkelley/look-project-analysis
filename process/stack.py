#!/usr/bin/env python3
import os
import json
import argparse
import warnings
from glob import glob
from itertools import groupby
from collections import defaultdict
from astropy.nddata import CCDData
import numpy as np
from scipy.cluster import hierarchy
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.table import vstack
from astropy.wcs import WCS, FITSFixedWarning
import astropy.units as u
from reproject import reproject_adaptive
from ccdproc import Combiner
from sbpy.data import Ephem, Names, TimeScaleWarning
from sbpy.data.names import TargetNameParseError
from quick_look import (
    assumed_gmr,
    locations,
    rho_labels,
    filters as filter_to_ps1,
    target_to_filename,
    setup_logger,
)

warnings.simplefilter("ignore", (UserWarning, FITSFixedWarning, TimeScaleWarning))


class NoGoodData(Exception):
    pass


def combine_photometry(phot, target, site, filter):
    i = (phot["merr5"] < 0.2) * np.isfinite(phot["m5"])

    # anything more than 3 pix away is probably bad
    if len(phot[i]) > 2:
        # if there are more than 2 detections, use the median as a reference
        cm = np.median(phot["dc"][i])
    else:
        # otherwise use the ephemeris position
        cm = 0
    i = i * (np.abs(phot["dc"].data - cm) < 3)

    if not any(i):
        raise NoGoodData()

    t = Time(phot["date"][i])
    jd = t.jd.mean()
    m = {}
    merr = {}
    for r in rho_labels:
        k = "{}".format(r)
        _m, sw = np.average(
            phot["m" + k][i], weights=phot["merr" + k][i] ** -2, returned=True
        )
        m["m" + k] = _m
        merr["merr" + k] = sw**-0.5

    row = {
        "target": target,
        "source": site,
        "date": Time(jd, format="jd").iso[:-4],
        "fracday": t.mjd.mean() % 1,
        "tmtp": phot["tmtp"].mean(),
        "ra": phot["ra"].mean(),
        "dec": phot["dec"].mean(),
        "dc": phot["dc"][i].mean(),
        "rh": phot["rh"][i].mean(),
        "delta": phot["delta"][i].mean(),
        "phase": phot["phase"][i].mean(),
        "pixel scale": phot["pixel scale"][i].mean(),
        "filter": filter,
        "catalog filter": phot["catalog filter"][0],
        "N exp": i.sum(),
        "exptime": phot["exptime"][i].sum(),
        "airmass": phot["airmass"][i].mean(),
        "seeing": phot["seeing"][i].mean(),
        "rho10k": phot["rho10k"][i].mean(),
        "meanzp": phot["zp"][i].mean(),
    }
    row.update(m)
    row.update(merr)

    return row


def color_and_outbursts(phot):
    pass


def grouper(row):
    return (row["target"], row["site"], row["filter"])


def stack_cluster(target, cluster, shape, site):
    dates = Time(cluster["mjd"], format="mjd")
    date = Time(dates.mjd.mean(), format="mjd")

    # create WCS objects in the comet's rest frame for each image and the stacked image
    opts = dict(epochs=dates, location=locations[site[:3]])
    try:
        target_type = Names.asteroid_or_comet(target)
    except TargetNameParseError:
        target_type = "asteroid"

    if target_type == "comet":
        id_type = "designation"
        opts["closest_apparition"] = True
        opts["no_fragments"] = True
    else:
        id_type = "smallbody"

    eph = Ephem.from_horizons(target, id_type, **opts)

    wcs = []
    for i, f in enumerate(cluster["file"]):
        w = WCS(fits.getheader(f, 1))
        w.wcs.crval[0] += (eph["ra"][0] - eph["ra"][i]).to_value("deg")
        w.wcs.crval[1] += (eph["dec"][0] - eph["dec"][i]).to_value("deg")
        wcs.append(w)

    wcs0 = WCS()
    wcs0.wcs.radesys = "ICRS"
    wcs0.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs0.wcs.crval = cluster["ra"][0], cluster["dec"][0]
    wcs0.wcs.crpix = shape // 2
    ps0 = u.Quantity(wcs[0].proj_plane_pixel_scales()).mean().to_value("arcsec")
    cd = wcs[0].wcs.cd * 0.39 / ps0

    # align with N up, E left
    # For E left:[0, 0] must be < 0
    # For N up: [1, 1] must be > 0
    if cd[0, 0] > 0:
        cd[0, 0] *= -1
    if cd[1, 1] < 0:
        cd[1, 1] *= -1
    wcs0.wcs.cd = cd

    wcs0.wcs.dateobs = date.isot
    wcs0.wcs.mjdobs = date.mjd

    # background subtract, photometric calibration, re-projection
    stack = np.ma.empty((len(cluster["file"]), shape[0], shape[1]))
    stack.mask = np.zeros(stack.shape, bool)
    scale = np.empty(len(stack))
    for i in range(len(cluster["file"])):
        im0 = fits.getdata(cluster["file"][i])
        im0 -= fits.getdata(
            f"backgrounds/{os.path.basename(cluster['file'][i][:-3])}.gz"
        )
        mzp = cluster["zp"][i] + cluster["color cor"][i] * assumed_gmr[target]
        scale[i] = 10 ** (-0.4 * (mzp - 25))
        im0 *= scale[i]

        # I'm not sure how best to avoid interpolation artifacts.  It may be because
        # the rotation is usually very small.  Using nearest-neighbor here looks good.
        # but bicubic is best for actual interpolations.
        # stack.data[i], cov = reproject_adaptive((im0, wcs[i]), wcs0, shape,
        #                                     order='nearest-neighbor')
        # ...
        # stack.data[i], cov = reproject_exact((im0, wcs[i]), wcs0, shape, parallel=True)
        stack.data[i], cov = reproject_adaptive((im0, wcs[i]), wcs0, shape, parallel=True, conserve_flux=True)
        stack.mask[i] = cov == 0

    stack = CCDData(stack, unit="adu")
    combiner = Combiner(stack)

    med = combiner.median_combine().data
    if len(stack.data) > 2:
        # median and average are the same for n = 2
        avg = combiner.average_combine().data
    else:
        avg = None

    return wcs0, avg, med, scale


def save_stack(wcs, avg, med, scale, target, site, cluster, filter, prefix):
    dates = Time(cluster["mjd"], format="mjd")
    date = Time(dates.mjd.mean(), format="mjd")

    hdu = wcs.to_fits()
    hdu[0].header["target"] = target
    hdu[0].header["site"] = site
    hdu[0].header["startmjd"] = (dates - cluster["exptime"] * u.s).min().mjd
    hdu[0].header["midmjd"] = date.mjd
    hdu[0].header["stopmjd"] = (dates + cluster["exptime"] * u.s).max().mjd
    hdu[0].header["filter"] = filter, "observation filter"
    hdu[0].header["calfilt"] = (
        filter_to_ps1[filter],
        "data is calibrated to this PS1 filter",
    )
    hdu[0].header["gain"] = scale.mean() ** -1, "mean eff. gain per image"
    hdu[0].header["ra"] = cluster["ra"].mean(), "mean observed RA, deg"
    hdu[0].header["dec"] = cluster["dec"].mean(), "mean observed Dec, deg"
    hdu[0].header["rh"] = cluster["rh"].mean()
    hdu[0].header["delta"] = cluster["delta"].mean()
    hdu[0].header["phase"] = cluster["phase"].mean()
    hdu[0].header["rho10k"] = cluster["rho10k"].mean()
    hdu[0].header["airmass"] = cluster["airmass"].mean()
    hdu[0].header["seeing"] = (
        cluster["seeing"].mean(),
        "mean approximate seeing, arcsec",
    )

    try:
        hdu[0].header["tmtp"] = cluster["tmtp"].mean(), "T-Tp"
    except ValueError:
        # in case tmtp is masked
        pass
    hdu[0].header["exptime"] = cluster["exptime"].sum()
    hdu[0].header["mzp"] = 25.0, "magnitude zero point"
    hdu[0].header["mzpunc"] = (
        np.sqrt(np.sum(cluster["zp err"] ** 2)) / np.sqrt(len(cluster)),
        "zeropoint uncertainty",
    )
    hdu[0].header["colorcor"] = (cluster["color cor"].mean(), "mean color correction")
    hdu[0].header["tgtcolor"] = assumed_gmr[target], "target g-r used for calibration"

    hdu[0].header["nimages"] = len(cluster), "number of images combined"
    hdu[0].header.add_comment(
        "Files combined: {}".format([os.path.basename(f) for f in cluster["file"]])
    )

    stack_files = []
    for im, label in zip((avg, med), ("avg", "med")):
        if im is None:
            continue
        hdu[0].header["combfunc"] = label, "combination function"
        hdu[0].data = im.astype(np.float32)
        stack_files.append(f"{prefix}_{label}.fits")
        hdu.writeto(stack_files[-1], overwrite=True)

    return stack_files


def update_binned_photometry(cluster, target, site, filter, prefix, binned, phot):
    try:
        row = combine_photometry(cluster, target, site, filter)
    except NoGoodData:
        return

    row["stack prefix"] = prefix
    binned[prefix] = row

    tab = Table(rows=list(binned.values()), meta=phot.meta)
    tab["fracday"].format = "{:.5f}"
    tab["tmtp"].format = "{:.3f}"
    tab["ra"].format = "{:.6f}"
    tab["dec"].format = "{:.6f}"
    tab["dc"].format = "{:.1f}"
    tab["rh"].format = "{:.3f}"
    tab["delta"].format = "{:.3f}"
    tab["phase"].format = "{:.3f}"
    tab["pixel scale"].format = "{:.3f}"
    tab["exptime"].format = "{:.0f}"
    tab["airmass"].format = "{:.3f}"
    tab["seeing"].format = "{:.2f}"
    tab["rho10k"].format = "{:.2f}"
    tab["meanzp"].format = "{:.2f}"
    # tab['ostat'].format = '{:.1f}'
    for r in rho_labels:
        tab["m{}".format(r)].format = "{:.3f}"
        tab["merr{}".format(r)].format = "{:.3f}"

    tab.write("phot-binned.txt", format="ascii.fixed_width_two_line", overwrite=True)


def load_photometry():
    phot = ascii.read("phot.txt")
    phot["mjd"] = Time(phot["date"]).mjd
    phot["site"] = [f.split("/")[2].split("-")[0] for f in phot["file"]]

    phot.sort("date")
    phot.sort(("target", "site", "filter"))
    return phot


def load_binned_photometry():
    """inverted dictionary to key by stack prefix"""
    binned = {}
    if os.path.exists("phot-binned.txt"):
        for row in ascii.read("phot-binned.txt"):
            binned[row["stack prefix"]] = {
                k: v for (k, v) in zip(row.keys(), row.values())
            }
    return binned


def main(args, logger):
    phot = load_photometry()
    targets = set(phot["target"]) if args.target is None else [args.target]
    filters = set(phot["filter"]) if args.filter is None else [args.filter]
    logger.debug("Targets: %s", ", ".join(targets))
    logger.debug("Filters: %s", ", ".join(filters))

    binned = load_binned_photometry()

    shape = np.array((int(args.size * 60 / 0.39), int(args.size * 60 / 0.39)))

    associations = {}
    if os.path.exists(args.cluster_file):
        with open(args.cluster_file, "r") as inf:
            associations = json.load(inf)

    existing_files = set(glob("stacks/*/*.fits"))

    n = 0
    for (target, site, filter), rows in groupby(phot, grouper):
        if (target not in targets) or (filter not in filters):
            continue

        logger.debug("%s %s %s", target, site, filter)

        # this vstack is too slow: group = vstack(list(group))
        i = np.array([row.index for row in rows])
        group = phot[i]

        target_fn = target_to_filename(target)
        if not os.path.exists(f"stacks/{target_fn}"):
            os.system(f"mkdir -p stacks/{target_fn}")

        if len(group) > 1:
            clusters = hierarchy.fclusterdata(
                group["mjd"][:, np.newaxis], args.cluster_dt / 24, criterion="distance"
            )
        else:
            clusters = np.array([1])

        for cluster in [group[clusters == c] for c in set(clusters)]:
            cluster.sort("date")
            dates = Time(cluster["mjd"], format="mjd")
            date = Time(dates.mjd.mean(), format="mjd")

            if args.date is not None and date.iso[:10] != args.date:
                continue

            date_fn = date.isot.replace("-", "").replace(":", "")[:13]
            prefix = f"stacks/{target_fn}/{target_fn}_{date_fn}_{cluster['rh'].mean():.3f}au_{filter}_look"

            logger.debug("checking %s", prefix)
            process_stack = (
                # user requested all files to be processed?
                args.force
                # missing both stack files?
                or not any(
                    f"{prefix}_{label}.fits" in existing_files
                    for label in ("avg", "med")
                )
                # any files from this cluster not associated with this prefix?
                or any(
                    not associations.get(f, [""])[0].startswith(prefix)
                    for f in cluster["file"]
                )
            )
            if process_stack:
                logger.info("process %s", prefix)

                # Were these files already in a stack, even one with a different
                # prefix than being considered here?  Delete the stack and
                # revise the associations.
                old_stacks = set()

                for f in cluster["file"]:
                    if f in associations:
                        old_stacks = old_stacks.union(associations[f])
                        del associations[f]

                for old_stack in old_stacks:
                    if os.path.exists(old_stack):
                        if args.process:
                            os.unlink(old_stack)
                            # also remove the temporally filtered image
                            try:
                                os.unlink(old_stack.replace(".fits", ".png"))
                            except FileNotFoundError:
                                pass
                        existing_files.discard(old_stack)
                    old_stack_prefix = os.path.basename(old_stack)[:-9]
                    if old_stack_prefix in binned:
                        del binned[old_stack]

                if args.process:
                    wcs, avg, med, scale = stack_cluster(target, cluster, shape, site)

                    stack_files = save_stack(
                        wcs, avg, med, scale, target, site, cluster, filter, prefix
                    )
                else:
                    # dummy update
                    stack_files = [f"{prefix}_med.fits"]

                for f in stack_files:
                    existing_files.add(f)

                # Update association file.
                for f in cluster["file"]:
                    associations[f] = stack_files
                if args.process:
                    with open(args.cluster_file, "w") as outf:
                        json.dump(associations, outf)

                n += 1

            # binned photometry should always be keyed by prefix
            if process_stack or prefix not in binned:
                if args.process:
                    update_binned_photometry(
                        cluster, target, site, filter, prefix, binned, phot
                    )

            if n >= args.n:
                return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-processing",
        dest="process",
        action="store_false",
        help="only identify files to be updated",
    )
    parser.add_argument("--target", help="only stack this target")
    parser.add_argument("--date", help="only stack on this UT date")
    parser.add_argument(
        "--cluster-dt",
        type=float,
        default=1,
        help="distance threshold for clustering (hours)",
    )
    parser.add_argument(
        "--cluster-file",
        default="stack-clusters.json",
        help="save cluster associations to this file",
    )
    parser.add_argument("--filter", help="only stack this filter")
    parser.add_argument(
        "--force", "-f", action="store_true", help="overwrite existing images"
    )
    parser.add_argument(
        "--size", "-s", default=7, type=int, help="image dimensions, arcmin"
    )
    parser.add_argument("-n", type=int, default=2**16, help="quit after `n` stacks")
    parser.add_argument(
        "-v", dest="level", action="store_const", const="DEBUG", help="verbose mode"
    )
    parser.add_argument(
        "-q", dest="level", action="store_const", const="WARNING", help="quiet mode"
    )
    args = parser.parse_args()

    logger = setup_logger(__file__, args.level)
    logger.info("Start")
    if not args.process:
        logger.info("No process mode enabled.  Stacks will not be generated.")

    try:
        main(args, logger)
    except:
        logger.exception("Unexpected exception occurred!")

    logger.info("Finish")
