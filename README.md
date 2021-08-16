# thinslice

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dunetpc v09_28_00 -q e20:prof

setup cmake v3_19_6

git clone https://gitlab.cern.ch/RooUnfold/RooUnfold.git

mkdir RooUnfold/build

cd RooUnfold/build

cmake ..

make -j4

cd ../..

git clone https://github.com/yangtj207/thinslice.git

cd thinslice

mkdir build

mkdir install

cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install

make install -j4

