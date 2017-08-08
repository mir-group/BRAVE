#!/bin/bash

MPI="mpirun"
VASP="vasp.5.4.4/bin/vasp_std"
PYTHON="python3"
POTPAW="potpaw_LDA.52"

cd 1_scf
cp $POTPAW/Si/POTCAR .
$MPI $VASP > silicon.vasp.out
cd ../2_bands
$PYTHON kpath.py
cp ../1_scf/POTCAR .
cp ../1_scf/POSCAR .
cp ../1_scf/CHGCAR .
$MPI $VASP > silicon.vasp.out
