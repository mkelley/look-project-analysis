#!/usr/bin/env python3
import json
import numpy as np
from astropy.io import ascii
from astropy.time import Time
import astropy.units as u
from mskpy import leading_num_key


def validate(value):
    if isinstance(value, np.int64):
        value = int(value)

    if not isinstance(value, str):
        return value if np.isfinite(value) else None

    return value


def serialize(tab):
    data = []
    for row in tab:
        data.append({k: validate(v) for k, v in zip(row.colnames, row.values())})
    return data


def json_format(tab):
    data = serialize(tab)
    return json.dumps(data)


phot = ascii.read("phot.txt")
with open("phot.json", "w") as outf:
    outf.write(json_format(phot))

binned = ascii.read("phot-binned.txt")
with open("phot-binned.json", "w") as outf:
    outf.write(json_format(binned))

targets = sorted(set(binned["target"]), key=leading_num_key)
rows = []
for target in targets:
    data = binned[binned["target"] == target]
    rows.append(
        {
            "target": target,
            "dates": (min(data["date"])[:10], max(data["date"])[:10]),
            "tmtp": (min(data["tmtp"]), max(data["tmtp"])),
            "rh": (min(data["rh"]), max(data["rh"])),
            "nobs": len(data),
            "nr": int((data["filter"] == "rp").sum()),
            "ng": int((data["filter"] == "gp").sum()),
        }
    )

recent = binned["date"] > (Time.now() - 1 * u.day).iso

summary = {
    "ntargets": len(rows),
    "nexposures": len(phot),
    "nobs": len(binned),
    "recent": serialize(binned[recent]),
    "targets": rows,
}

with open("summary.json", "w") as outf:
    json.dump(summary, outf)
