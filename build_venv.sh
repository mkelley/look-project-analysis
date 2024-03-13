#!/bin/bash -e

if [ ! -e .venv ]; then
    python3 -m venv .venv --prompt=LOOK
fi

source .venv/bin/activate

pip install -U pip
pip install -U setuptools wheel
pip install "git+https://github.com/mkelley/mskpy.git#egg=mskpy"
pip install "git+https://github.com/dfm/casjobs"
pip install -e .
cd /oort/msk/lco/look
wget https://raw.githubusercontent.com/mkelley/lco-sync/main/sync.py
