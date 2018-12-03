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

#
# 1) Download the SMASH code
#    Currently it only works if one has an access to smash repository,
#    which is not public. In future the repository has to be replaced to
#    a public one.
git clone hyihp-repos@fias.uni-frankfurt.de:scm/smash.git smash/smash_code

cd smash/smash_code

#
# 2) Compile SMASH
#
mkdir build
cd build
cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
number_of_cores=`nproc --all`
echo "Compiling SMASH using ${number_of_cores} cores."
make -j${number_of_cores}