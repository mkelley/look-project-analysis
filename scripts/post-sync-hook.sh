#!/bin/bash -eu
LOGROTATE=/usr/sbin/logrotate
source /home/msk/lco/.venv/bin/activate

cd /oort/msk/lco/look
$LOGROTATE look-project-analysis/scripts/logrotate.config -s logrotate.state

python3 phot.py -q &&\
python3 stack.py -q &&\
python3 color-and-ostat.py -q &&\
python3 phot-to-json.py &&\
python3 temporal-filter.py -q

cp -f \
   colors.txt avg-colors.txt\
   phot.txt phot-binned.txt\
   phot.json phot-binned.json\
   summary.json stack-clusters.json\
   color-histogram.png\
   /home/msk/public_html/look/

cp -f \
   colors.txt avg-colors.txt\
   phot.txt phot-binned.txt\
   phot.json phot-binned.json\
   summary.json stack-clusters.json\
   color-histogram.png\
   /home/msk/public_html/look-browser/data

