"""This module defines common parameters and functions."""

import numpy as np

# These are the "2014 CODATA recommended values" taken from
# "The NIST Reference on Constants, Units, and Uncertainty"
# http://physics.nist.gov/cuu/
#
BOHR = 0.52917721067  # Bohr radius, in Angstrom
RYDBERG = 13.605693009  # Rydberg constant times hc, in eV
HARTREE = 27.21138602  # Hartree energy, in eV
PLANCK = 4.135667662e-15  # Planck constant, in eV s
LIGHT = 2.99792458e8  # speed of light in vacuum, m / s
CHARGE = 1.6021766208e-19  # elementary charge, in Coulomb
MASS = 9.10938356e-31  # electron mass, in kg
BOLTZMANN = 8.6173303e-5  # Boltzmann constant, in eV / K

EPS12 = 1.0e-12
INF12 = 1.0e12

_ascale = {
        'bohr': 1.0 / BOHR,
        'angstrom': 1.0,
        'nm': 0.1}

_a3scale = {
        'uc': 0.0,
        'bohr3': _ascale['bohr'] ** 3,
        'angstrom3': _ascale['angstrom'] ** 3,
        'nm3': _ascale['nm'] ** 3}

_kscale = {
        'cartesian': 0.0,
        'crystal': 0.0}

_escale = {
        'ev': 1.0,
        'rydberg': 1.0 / RYDBERG,
        'hartree': 1.0 / HARTREE,
        'thz': 1.0e-12 / PLANCK,
        'cm-1': 1.0e-2 / PLANCK / LIGHT}

def _int_pts(value, param):
    dummy = np.abs(value - param)
    idx1 = dummy.argmin()
    dummy[idx1] = INF12
    idx2 = dummy.argmin()

    if (value[idx1] - param) * (value[idx2] - param) < 0.0:
        weight1 = abs(value[idx2] - param) / (
                abs(value[idx1] - param) + abs(value[idx2] - param))
        weight2 = abs(value[idx1] - param) / (
                abs(value[idx1] - param) + abs(value[idx2] - param))
    else:
        weight1 = 1.0
        weight2 = 0.0

    return idx1, idx2, weight1, weight2

