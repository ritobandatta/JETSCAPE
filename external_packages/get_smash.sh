#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# 1) Download the SMASH code
git clone --depth=1 https://github.com/smash-transport/smash.git --branch SMASH-3.0 smash/smash_code

# 2) Compile SMASH
cd smash/smash_code
mkdir build
cd build
#cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
cmake .. -DPythia_CONFIG_EXECUTABLE=/wsu/home/hl/hl97/hl9735/pythia8309/bin/pythia8-config -DCMAKE_PREFIX_PATH=/wsu/home/hl/hl97/hl9735/eigen-3.3.9  -DGSL_ROOT_DIR=/wsu/home/hl/hl97/hl9735/gsl
num_cores=${1:-1}
echo "Compiling SMASH using ${num_cores} cores."
make -j${num_cores} smash_shared

