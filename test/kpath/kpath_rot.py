# This script reads the k-points from files fcc.full.in and bcc.full.in,
# rotates them to the lattice vectors used by Quantum ESPRESSO (ibrav = 2 & 3),
# and writes them to files fcc_rot.full.in and bcc_rot.full.in.

import numpy as np
import brave

kpt = brave.Kpoint()
kpt.read('internal', ['fcc.full.in'])
del kpt.bvec
avec_rot = np.array([[-0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.5, 0.0]], float)
mtx = np.dot(avec_rot, np.linalg.inv(kpt.avec))
kpath_rot = np.transpose(np.dot(mtx, np.transpose(kpt.kpath)))
kpt.avec = avec_rot
kpt.kpath = kpath_rot
kpt.write('internal', ['fcc_rot.full.in'])

kpt = brave.Kpoint()
kpt.read('internal', ['bcc.full.in'])
del kpt.bvec
avec_rot = np.array([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]], float)
mtx = np.dot(avec_rot, np.linalg.inv(kpt.avec))
kpath_rot = np.transpose(np.dot(mtx, np.transpose(kpt.kpath)))
kpt.avec = avec_rot
kpt.kpath = kpath_rot
kpt.write('internal', ['bcc_rot.full.in'])
