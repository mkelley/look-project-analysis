#!/bin/bash -eu
LOGROTATE=/usr/sbin/logrotate
source /home/msk/lco/.venv/bin/activate
source /oort/msk/lco/look/google-cloud-sdk/path.bash.inc

cd /oort/msk/lco/look
$LOGROTATE look-project-analysis/scripts/logrotate.config -s logrotate.state

python3 file-summary.py &&\
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

zip phot.json.zip phot.json &&\
zip phot-binned.json.zip phot-binned.json &&\
zip summary.json.zip summary.json &&\
zip stack-clusters.json.zip stack-clusters.json &&\
zip phot.txt.zip phot.txt &&\
zip phot-binned.txt.zip phot-binned.txt &&\
gsutil cp phot.json.zip gs://lco-outbursting-objects-0514.appspot.com/data/phot.json.zip &&\
gsutil cp phot-binned.json.zip gs://lco-outbursting-objects-0514.appspot.com/data/phot-binned.json.zip &&\
gsutil cp summary.json.zip gs://lco-outbursting-objects-0514.appspot.com/data/summary.json.zip &&\
gsutil cp stack-clusters.json.zip gs://lco-outbursting-objects-0514.appspot.com/data/stack-clusters.json.zip &&\
gsutil cp phot.txt.zip gs://lco-outbursting-objects-0514.appspot.com/data/phot.txt.zip &&\
gsutil cp phot-binned.txt.zip gs://lco-outbursting-objects-0514.appspot.com/data/phot-binned.txt.zip
