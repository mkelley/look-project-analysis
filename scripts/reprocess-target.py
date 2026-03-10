
#!/usr/bin/env python3
import os
import argparse
import subprocess
from glob import glob
from astropy.io import ascii
from astropy.table import Table
import quick_look as ql

parser = argparse.ArgumentParser()
parser.add_argument("target", help="reprocess this target name and its aliases")
args = parser.parse_args()

# get all target aliases
target_names = [args.target] + [
    k for k, v in ql.rename_target.items() if v == args.target
]

print("Processing targets:")
print("  - ", "  - \n".join(target_names), sep="")

# get all files labeled with these targets
summary = ascii.read("file-summary.txt")

reprocess = set()
for row in summary:
    if row["target"] in target_names:
        reprocess.add(row["filename"])

if not any(reprocess):
    print("No files to reprocess.")
    exit()

# remove from the skip list
if os.path.exists("phot-skip.list"):
    subprocess.check_call(["cp", "phot-skip.list", "phot-skip.list~"])
    skip = ascii.read("phot-skip.list")

    # skip list is just basename sans extension
    _reprocess = {os.path.basename(fn)[:-8] for fn in reprocess}
    keep = []
    for row in skip:
        keep.append(row["file"] not in _reprocess)

    if all(keep):
        print("phot-skip.list unchanged")
    else:
        count = len(skip) - sum(keep)
        skip = skip[keep]
        skip.write("phot-skip.list", format="ascii.csv", overwrite=True)
        print(f"{count} rows removed from phot-skip.list")
else:
    print("phot-skip.list does not exist")

# remove from the photometry results
if os.path.exists("phot.txt"):
    subprocess.check_call(["cp", "phot.txt", "phot.txt~"])
    phot = ql.set_phot_table_formats(ascii.read("phot.txt"))

    keep = []
    for row in phot:
        keep.append(row["file"] not in reprocess)

    if all(keep):
        print("phot.txt unchanged")
    else:
        count = len(phot) - sum(keep)
        phot = phot[keep]
        phot.write("phot.txt", format="ascii.fixed_width_two_line", overwrite=True)
        print(f"{count} rows removed from phot.txt")
else:
    print("phot.txt does not exist")

# remove all backgrounds
count = 0
for fn in reprocess:
    bgf = f"backgrounds/{os.path.basename(fn)}"
    if os.path.exists(bgf):
        os.unlink(bgf)
        count += 1
print(f"Removed {count} background files.")

# remove all stacks
target_fn = ql.target_to_filename(args.target)
files = glob(f"stacks/{target_fn}/{target_fn}-*.fits")
for fn in files:
    os.unlink(fn)
print(f"Removed {len(files)} stack files.")

# remove all processed data
files = glob(f"processed/{target_fn}/{target_fn}-*.fits")
for fn in files:
    os.unlink(fn)
print(f"Removed {len(files)} processed files.")

files = glob(f"plots/{target_fn}/{target_fn}-*.jpg")
for fn in files:
    os.unlink(fn)
print(f"Removed {len(files)} plots.")

print("Ready to reprocess (e.g., with post-sync-hook.sh)")
