#!/bin/bash

echo "Building the project..."

if [ ! -d "./build" ]; then
    mkdir build
fi

cd build

if cmake ../ -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZER_ADDRESS=True $@; then
echo "calling make...";
else 
echo "CMake failed!";
cd ..
exit 1 
fi

num_procs=$(nproc)
echo "building using ${num_procs} workers"
if make -j${num_procs}; then
echo "build ready";
else 
echo "build failed!";
cd ..
exit 2
fi

cd ..

