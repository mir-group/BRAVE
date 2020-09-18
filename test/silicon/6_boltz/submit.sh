#!/bin/bash

PYTHON="python3"
BOLTZTRAP="boltztrap-1.2.5/src/BoltzTraP"

cp ../1_qe/3_ph/silicon.epa.e .
$PYTHON qe2boltz.py > qe2boltz.out
$BOLTZTRAP silicon.def > silicon.boltztrap.out
