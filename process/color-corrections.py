"""
Fit a color correction to all data.

This script assumes data already processed by phot.py (i.e., the RefCat2
database is already populated with relevant sources).
"""

import os
import numpy as np
from glob import glob
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from calviacat import RefCat2
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
)

cctabfn = "color-corrections.csv"
if os.path.exists(cctabfn):
    tab = ascii.read(cctabfn)
    rows = [dict(zip(tab.colnames, row)) for row in tab]
else:
    rows = []

skip = set(ascii.read("phot-skip.list")["file"])
skip |= set([row["file"] for row in rows])

rc2 = RefCat2("cat.db", min_matches=10)

files = [fn for fn in glob("e91/*/*fz")
         if os.path.basename(fn)[:-8] not in skip]

for fn in files:
    basename = os.path.basename(fn)[:-8]

    with fits.open(fn) as hdul:
        h = hdul["sci"].header
        phot = Table(hdul["cat"].data)

    phot = phot[phot["flux"] > 0]
    lco = SkyCoord(phot["ra"], phot["dec"], unit="deg")
    try:
        objids, distances = rc2.xmatch(lco)
    except TypeError:
        row = {
            "file": basename,
            "n": 0
        }

    m_inst = -2.5 * np.log10(phot["flux"])
    m_err = phot["fluxerr"] / phot["flux"] * 1.0857
    catfilt = filters.get(h["FILTER"])
 
    try:
        zp, C, zp_unc, m, gmr, gmi = rc2.cal_color(
            objids, m_inst, catfilt, "g-r", mlim=[13, 18], gmi_lim=[0.2, 3.0]
        )
        n_cal = (~m.mask).sum()
        row = {
            "file": basename,
            "siteid": h["SITEID"],
            "encid": h["ENCID"],
            "telid": h["TELID"],
            "airmass": h["AIRMASS"],
            "relhum": h["WMSHUMID"],
            "obsfilt": h["FILTER"],
            "catfilt": catfilt,
            "zp": zp,
            "C": C,
            "zp_unc": zp_unc,
            "n": (~m.mask).sum(),
        }
    except Exception:
        row = {
            "file": basename,
            "n": 0
        }
    rows.append(row)

cctab = Table(rows)
os.system(f"cp --backup=numbered {cctabfn} {cctabfn}")
cctab.write(cctabfn, overwrite=True)