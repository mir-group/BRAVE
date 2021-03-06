#!/bin/bash

WANNIER90="wannier90-2.1.0/wannier90.x"
POSTW90="wannier90-2.1.0/postw90.x"

$WANNIER90 -pp silicon
mv silicon.wout silicon.wout_1
cp ../1_qe/5_wan/silicon.amn .
cp ../1_qe/5_wan/silicon.eig .
cp ../1_qe/5_wan/silicon.mmn .
$WANNIER90 silicon
mv silicon.wout silicon.wout_2
$POSTW90 silicon
