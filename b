#!/bin/bash -v
mkdir -p build && cd build && cmake .. && make -I/usr/include/coin
