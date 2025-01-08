#!/bin/bash
git clone --recursive https://github.com/mcale6/dr_sasa_python.git 
cd dr_sasa_python

python3.10 -m venv dr_sasa_venv
source dr_sasa_venv/bin/activate

pip install -e .

cmake -B build && cmake --build build && cmake --install build

export PYTHONPATH="$PWD/dr_sasa_venv/lib/python3.10/site-packages:$PWD/build:$PYTHONPATH"

pytest -v