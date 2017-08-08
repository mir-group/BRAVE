# This script reads the k-points from files fcc.full.in and bcc.full.in,
# rotates them to the lattice vectors used by Quantum ESPRESSO (ibrav = 2 & 3),
# and writes them to files fcc_rot.full.in and bcc_rot.full.in.

import numpy
import brave

kpt = brave.Kpoint()
kpt.read('internal', ['fcc.full.in'])
del kpt.bvec
avec_rot = numpy.array([[-0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.5, 0.0]], float)
mtx = numpy.dot(avec_rot, numpy.linalg.inv(kpt.avec))
kpath_rot = numpy.transpose(numpy.dot(mtx, numpy.transpose(kpt.kpath)))
kpt.avec = avec_rot
kpt.kpath = kpath_rot
kpt.write('internal', ['fcc_rot.full.in'])

kpt = brave.Kpoint()
kpt.read('internal', ['bcc.full.in'])
del kpt.bvec
avec_rot = numpy.array([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]], float)
mtx = numpy.dot(avec_rot, numpy.linalg.inv(kpt.avec))
kpath_rot = numpy.transpose(numpy.dot(mtx, numpy.transpose(kpt.kpath)))
kpt.avec = avec_rot
kpt.kpath = kpath_rot
kpt.write('internal', ['bcc_rot.full.in'])
