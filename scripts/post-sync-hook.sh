#!/bin/bash -eu
source /home/msk/lco/.venv/bin/activate
cd /oort/msk/lco/look
python3 phot.py
for target in 2005qn173 29p c2014un271 382p 57p 97p 67p 191p; do
    PYTHONPATH=. python3 look-project-analysis/figures/plot-${target}.py
done
#python3 summarize-for-web.py
python3 stack.py
cp -f avg-colors.txt\
   c2014un271-look-phot-binned.txt\
   colors.txt\
   phot-binned-date-sort.txt\
   phot.txt phot-binned.txt\
   /home/msk/public_html/look/
rsync 29p-look.png\
      2005qn173-look.png\
      141p-look.png\
      156p-look.png\
      7p-look-periodigram.png\
      7p-look-phased.png\
      7p-look.png\
      c2014un271-look.png\
      c2020r4-look.png\
      382p-look.png\
      57p-look.png\
      97p-look.png\
      191p-look.png\
      67p-look.png\
      7p-look-20210710-outburst.gif\
      /home/msk/public_html/look/plots/


