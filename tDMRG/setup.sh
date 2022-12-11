#!/bin/bash

mkdir output
mkdir output/open
mkdir output/stop
wget 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz'
tar -xzf eigen-3.4.0.tar.gz
mv eigen-3.4.0 ../Eigen
rm eigen-3.4.0.tar.gz
