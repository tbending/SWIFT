#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f advection_1D.hdf5 ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift --hydro --radiation --stars --feedback --external-gravity rt_advection1D.yml

python3 plotSolution.py
