# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin
#export ROOTSYS=$HOME/root
#export PATH=$PATH:$ROOTSYS/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
#export PATH="~/cmake-3.4.1-Linux-x86_64/bin:$PATH"
#module load boost/1.67.0
#module load gsl/2.5
#module load hdf5/1.10.2 
#module load root/6.28.10 
module load cmake/3.21.1
#module load pythia/8.303
#module load openmpi3/3.1.0
#module load fftw/3.3.7
#module load eigen/3.3.7 
module load python/3.10

# Set up Pythia 8.309
export PYTHIA8_ROOT_DIR=/wsu/home/hl/hl97/hl9735/pythia8309
export PYTHIA8_DIR=${PYTHIA8_ROOT_DIR}
export PYTHIA8_INCLUDE_DIR=${PYTHIA8_DIR}/include
export PYTHIA8_LIB=${PYTHIA8_DIR}/lib
export PYTHIA8_XML=${PYTHIA8_DIR}/xmldoc

# Add Pythia binaries to PATH
export PATH=${PYTHIA8_DIR}/bin:${PATH}

# Add Pythia libraries to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${PYTHIA8_LIB}:${LD_LIBRARY_PATH}

# Verify the correct pythia8-config is found
#echo "Using pythia8-config at $(which pythia8-config)"

export PATH=/wsu/home/hl/hl97/hl9735/pythia8309/bin:$PATH
export LD_LIBRARY_PATH==/wsu/home/hl/hl97/hl9735/pythia8309/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=/wsu/home/hl/hl97/hl9735/pythia8309:$CMAKE_PREFIX_PATH


source /wsu/home/hl/hl97/hl9735/root_dir/bin/thisroot.sh

export LD_LIBRARY_PATH=/wsu/home/hl/hl97/hl9735/gsl/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/wsu/home/hl/hl97/hl9735/gsl/include/gsl:$C_INCLUDE_PATH
export GSL_HOME=/wsu/home/hl/hl97/hl9735/gsl

#export SMASH_DIR=/wsu/home/hl/hl97/hl9735/JETSCAPE/external_packages/smash/smash_code


export BOOST_ROOT=/wsu/home/hl/hl97/hl9735/boost_
export LD_LIBRARY_PATH=$BOOST_ROOT/lib:$LD_LIBRARY_PATH
export PATH=$BOOST_ROOT/bin:$PATH

export PATH
