#!/bin/bash

mkdir build
cp ./src/python_wrapper.py ./build
cd build
cmake -DCMAKE_Fortran_COMPILER=/usr/local/bin/gfortran ../
make clean; make
