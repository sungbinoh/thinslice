# thinslice

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dunetpc v09_25_00 -q e20:prof

setup cmake v3_19_6

git clone -b 2.1 https://gitlab.cern.ch/RooUnfold/RooUnfold.git

mkdir RooUnfold/build

cd RooUnfold/build

cmake ..

make -j4

cd ../..

git clone git@github.com:fnal-dunephysics/thinslice.git

cd thinslice

mkdir build

mkdir install

cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install

make install -j4

