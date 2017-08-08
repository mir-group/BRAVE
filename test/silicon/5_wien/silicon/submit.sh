#!/bin/bash

PYTHON="python3"

init_lapw -b -red 0 -vxc 5 -ecut -6.0 -rkmax 9.0 -numk 512 > silicon.outputinit
run_lapw -p
$PYTHON kpath.py
x_lapw lapw1 -band -p > silicon.outputlapw1
cat silicon.output1_? > silicon.output1
cat silicon.output1_?? >> silicon.output1
