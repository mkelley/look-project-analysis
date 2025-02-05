import os
import sys
from itertools import groupby
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from astropy.visualization import hist
from astropy.stats import sigma_clip
from mskpy import between

try:
    cc = ascii.read("color-corrections.csv")
except:
    print("Requires color-corrections.csv from color-corrections.py")
    sys.exit()

cc = cc[cc["n"] > 0]
cc.sort(["telid", "obsfilt", "catfilt"])
cc = [dict(zip(cc.colnames, row)) for row in cc]

summary = []
for group, rows in groupby(cc, lambda row: "-".join([row["telid"], row["obsfilt"], row["catfilt"]])):
    data = Table(list(rows))

    if len(data) < 100:
        continue

    i = between(data["C"], [-0.2, 0.3])

    plt.clf()
    hist(data["C"][i], 30)
    plt.setp(plt.gca(), ylabel="$N$", xlabel="$C$ (mag/mag)")
    plt.tight_layout(pad=0.2)
    plt.savefig(f"color-corrections-{group}.png", dpi=300)

    sc = sigma_clip(data["C"])
    summary.append({
        "group": group,
        "N": np.sum(~sc.mask),
        "mean": np.mean(sc),
        "median": np.median(sc),
        "stdev": np.std(sc),
    })

tab = Table(summary)
fn = "color-corrections-summary.txt"
#os.system(f"cp --backup=numbered {fn} {fn}")
tab.write(fn, format="ascii.fixed_width_two_line", overwrite=True)


