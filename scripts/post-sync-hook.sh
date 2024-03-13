#!/bin/bash -eu
LOGROTATE=/usr/sbin/logrotate
source /oort/msk/lco/look/look-project-analysis/.venv/bin/activate
source /oort/msk/lco/look/google-cloud-sdk/path.bash.inc

cd /oort/msk/lco/look
$LOGROTATE look-project-analysis/scripts/logrotate.config -s logrotate.state

file-summary.py &&\
phot.py -q &&\
stack.py -q &&\
color-and-ostat.py -q &&\
phot-to-json.py &&\
temporal-filter.py -q

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
