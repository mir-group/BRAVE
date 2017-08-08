This is the BRAVE package.

BRAVE stands for Bloch Representation Analysis and Visualization Environment.
It parses the output of various electronic structure codes, generates the input
files for subsequent calculations, and helps to analyse and plot the results.

For details please refer to the documentation on individual modules:
>>> help(brave.Energy)
>>> help(brave.Transport)
>>> help(brave.Plot)

List of the supported file formats and the corresponding codes:

fileformat       filename             executable       code       access
----------       --------             ----------       ----       ------
'internal'       'prefix.brave'         --               --       rw
'pw-in'          'prefix.in'          pw.x             QE [1]     w
'pw-out'         'prefix.out'         pw.x             QE [1]     r
'bands-out'      'bands.out'          bands.x          QE [1]     r
'matdyn-out'     'matdyn.modes'       matdyn.x         QE [1]     r
'matdyn-dos'     'prefix.vdos'        matdyn.x         QE [1]     r
'epa-out'        'prefix.epa'         epa.x            QE [1]     r
'inteqp-out'     'bandstructure.dat'  inteqp.flavor.x  BGW [2]    r
'sigma-out'      'sigma_hp.log'       sigma.flavor.x   BGW [2]    r
'wannier-in'     'seedname.win'       wannier90.x      Wannier90  rw
'wannier-out'    'seedname_band.dat'  wannier90.x      Wannier90  r
'vasp-kpt'       'KPOINTS'            vasp             VASP       w
'vasp-out'       'OUTCAR'             vasp             VASP       r
'lapw-kpt'       'case.klist_band'    lapw1            WIEN2k     w
'lapw-out'       'case.output1'       lapw1            WIEN2k     r
'boltztrap-in'   [3]                  BoltzTraP        BoltzTraP  w
'boltztrap-out'  [4]                  BoltzTraP        BoltzTraP  r
'boltztrap-dos'  [5]                  BoltzTraP        BoltzTraP  r

------------------------------------------------------------------------
[1] Quantum ESPRESSO
[2] BerkeleyGW
[3] 'prefix.def', 'prefix.intrans', 'prefix.struct', 'prefix.energy[so]'
[4] 'prefix.intrans', 'prefix.trace'
[5] 'prefix.intrans', 'prefix.transdos'

In the case of a spin-polarized system, some executables
store the two spin components in separate files. These are
listed below, following the naming convention of WIEN2k:

fileformat     filename
----------     --------
'bands-out'    'bands.out[up|dn]'
'wannier-out'  'seedname_band.dat[up|dn]'
'lapw-out'     'case.output1[up|dn]'

Notes on different codes:

Quantum ESPRESSO: pw.x must be run with verbosity = 'high'
to force writing symmetry operations and k-points to file
'prefix.out'.

Wannier90: Property kpoint of class Kpoint is not available
from the input or output files of wannier90.x. It can be
constructed from property kpath read from file 'seedname.win'
using methods calc_kpath and calc_kpoint.

VASP: Symmetry operations are not available from the input
or output files of vasp. They can be constructed from file
'POSCAR' using Atomic Simulation Environment and spglib as
done in BoltzTraP VASP interface. This is not implemented.

WIEN2k: Symmetry operations in cartesian coordinates are
available from file 'case.struct'. They can be read and
converted to crystal coordinates. This is not implemented.
File case.output1 can be merged together by hand:
    $ cat case.output1_? > case.output1
    $ cat case.output1_?? >> case.output1
    ... (etc. depending on how many files you have)
or generated automatically using spaghetti:
    $ x_lapw spaghetti [-p]

