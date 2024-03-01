import os
import json
from collections import defaultdict

# invert file -> cluster associations, but only track "median" files.
with open('stack-clusters.json', 'r') as inf:
    associations = json.load(inf)
stacks = defaultdict(set)
for fn, stack_files in associations.items():
    for stack in stack_files:
        if stack.endswith('med.fits'):
            stacks[stack].add(fn)

failed = []
passed = 0
for stack, sources in stacks.items():
    if not os.path.exists(stack):
        failed.append(stack)
    else:
        passed += 1
        for source in sources:
            if not os.path.exists(source):
                failed.append(source)
            else:
                passed += 1

print(f"{passed} files are present.")
print(f"{len(failed)} files are not present:")
print("\n".join(failed))

