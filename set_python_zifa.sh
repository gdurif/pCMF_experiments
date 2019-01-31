#!/bin/bash

#### set virtualenv and install ZIFA

## PROJECT ROOT DIR
PROJDIR=$(git rev-parse --show-toplevel)

## SET ENVIRONMENT
source $PROJDIR/set_env.sh

## virtual env
virtualenv -p $(which python3) ${PROJDIR}/.pyenv/zifa
source ${PROJDIR}/.pyenv/zifa/bin/activate

## required package
pip install numpy scipy matplotlib scikit-learn

## ZIFA
cd ${PROJDIR}/third_party_src/ZIFA
python setup.py install

## quit
cd ${PROJDIR}
deactivate
