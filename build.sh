export HPCKit_VERSION=$1
mkdir build
cd build
cmake ..
make -j
