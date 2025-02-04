#!/usr/bin/env python3
import os
from glob import glob
from astropy.io import ascii, fits
from astropy.time import Time
from astropy.table import Table


def row_to_dict(row):
    return {col: row[col] for col in row.colnames}


if os.path.exists("file-summary.txt"):
    tab = ascii.read("file-summary.txt")
    files = {row["filename"]: row_to_dict(row) for row in tab}
else:
    files = {}


for f in glob("e91/202*/*.fits.fz"):
    if f in files:
        continue

    try:
        h = fits.getheader(f, extname="sci")
    except OSError:
        row = {"filename": f}
        continue

    row = {
        "filename": f,
        "file time": Time(os.stat(f).st_mtime, format="unix").iso,
        "target": h["OBJECT"],
        "obs date": h["DATE-OBS"],
        "filter": h["FILTER"],
        "exp time": h["EXPTIME"],
        "userid": h["USERID"],
        "propid": h["PROPID"],
    }

    files[f] = row

tab = Table(sorted(files.values(), key=lambda row: row["obs date"]))
tab.sort("obs date")
tab.write("file-summary.txt", format="ascii.fixed_width_two_line", overwrite=True)
