#!/bin/bash

# epa.x dies with I/O errors when ph.x is run with buffered I/O
unset FORT_BUFFERED
unset FORT_BLOCKSIZE
unset FORT_BUFFERCOUNT

NPOOL=20
MPI="mpirun"
PW="q-e-qe-6.6/bin/pw.x"
BANDS="q-e-qe-6.6/bin/bands.x"
PH="q-e-qe-6.6/bin/ph.x"
Q2R="q-e-qe-6.6/bin/q2r.x"
MATDYN="q-e-qe-6.6/bin/matdyn.x"
EPA="q-e-qe-6.6/bin/epa.x"
PW2BGW="q-e-qe-6.6/bin/pw2bgw.x"
PW2WANNIER90="q-e-qe-6.6/bin/pw2wannier90.x"
PYTHON="python3"

cd 1_scf
$MPI $PW -npool $NPOOL < 1.pw.in > 1.pw.out
cd ../2_bands
cp -r ../1_scf/silicon.save .
$MPI $PW -npool $NPOOL < 1.pw.in > 1.pw.out
$MPI $BANDS -npool $NPOOL < 2.bands.in > 2.bands.out
rm silicon.wfc*
cd ../3_ph
cp -r ../1_scf/silicon.save .
$MPI $PH -npool $NPOOL < 1.ph.in > 1.ph.out
$MPI $PH -npool $NPOOL < 2.ph.in > 2.ph.out
cp silicon.dyn0 silicon.dyn0.xml
$MPI $Q2R < 3.q2r.in > 3.q2r.out
$MPI $MATDYN < 4.matdyn.in > 4.matdyn.out
$MPI $MATDYN < 5.matdyn.in > 5.matdyn.out
$EPA < 6.epa.in > 6.epa.out
cd ../4_bgw
$PYTHON kpath.py
cp -r ../1_scf/silicon.save .
$MPI $PW -npool $NPOOL < 1.pw.in > 1.pw.out
$MPI $PW2BGW -npool $NPOOL < 2.pw2bgw.in > 2.pw2bgw.out
$MPI $PW -npool $NPOOL < 3.pw.in > 3.pw.out
$MPI $PW2BGW -npool $NPOOL < 4.pw2bgw.in > 4.pw2bgw.out
$MPI $PW -npool $NPOOL < 5.pw.in > 5.pw.out
$MPI $PW2BGW -npool $NPOOL < 6.pw2bgw.in > 6.pw2bgw.out
rm silicon.wfc*
cd ../5_wan
cp -r ../1_scf/silicon.save .
cp ../../3_wan/silicon.nnkp .
$MPI $PW -npool $NPOOL < 1.pw.in > 1.pw.out
$MPI $PW2WANNIER90 < 2.pw2wannier90.in > 2.pw2wannier90.out
rm silicon.wfc*
cd ../6_boltz
cp -r ../1_scf/silicon.save .
$MPI $PW -npool $NPOOL < 1.pw.in > 1.pw.out
rm silicon.wfc*
