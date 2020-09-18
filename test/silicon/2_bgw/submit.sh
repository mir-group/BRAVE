#!/bin/bash

MPI="mpirun"
EPSILON="BerkeleyGW-2.1/bin/epsilon.real.x"
SIGMA="BerkeleyGW-2.1/bin/sigma.real.x"
INTEQP="BerkeleyGW-2.1/bin/inteqp.real.x"

cd 1_epsilon
cp ../../1_qe/4_bgw/WFN .
cp ../../1_qe/4_bgw/WFNq .
$MPI $EPSILON > epsilon.out
cd ../2_sigma
cp ../../1_qe/4_bgw/vxc.dat .
cp ../../1_qe/4_bgw/RHO .
cp ../../1_qe/4_bgw/WFN WFN_inner
cp ../1_epsilon/eps0mat .
cp ../1_epsilon/epsmat .
$MPI $SIGMA > sigma.out
cd ../3_inteqp
cp ../2_sigma/eqp1.dat eqp_co.dat
cp ../../1_qe/4_bgw/WFN WFN_co
cp ../../1_qe/4_bgw/WFN_fi .
$MPI $INTEQP > inteqp.out
