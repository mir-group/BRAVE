## Bloch Representation Analysis and Visualization Environment

BRAVE parses the output of different electronic structure codes, generates the input files for subsequent calculations, and helps to analyse and plot the results.

#### Supported codes and the corresponding arguments used by BRAVE

| Code                                                  | Executable          | Access | fileformat      | filenames                                                |
|-------------------------------------------------------|---------------------|--------|-----------------|----------------------------------------------------------|
|                                                       |                     | rw     | 'internal'      | 'prefix.brave'                                           |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **pw.x**            | w      | 'pw-in'         | 'prefix.in'                                              |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **pw.x**            | r      | 'pw-out'        | 'prefix.out'                                             |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **bands.x**         | r      | 'bands-out'     | 'bands.out'                                              |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **matdyn.x**        | r      | 'matdyn-out'    | 'matdyn.modes'                                           |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **matdyn.x**        | r      | 'matdyn-dos'    | 'prefix.vdos'                                            |
| [Quantum Espresso](https://www.quantum-espresso.org/) | **epa.x**           | r      | 'epa-out'       | 'epa.dat'                                                |
| [BerkeleyGW](https://www.berkeleygw.org/)             | **inteqp.flavor.x** | r      | 'inteqp-out'    | 'bandstructure.dat'                                      |
| [BerkeleyGW](https://www.berkeleygw.org/)             | **sigma.flavor.x**  | r      | 'sigma-out'     | 'sigma_hp.log'                                           |
| [Wannier90](https://www.wannier.org/)                 | **wannier90.x**     | rw     | 'wannier-in'    | 'seedname.win'                                           |
| [Wannier90](https://www.wannier.org/)                 | **wannier90.x**     | r      | 'wannier-out'   | 'seedname_band.dat'                                      |
| [VASP](https://www.vasp.at/)                          | **vasp**            | w      | 'vasp-kpt'      | 'KPOINTS'                                                |
| [VASP](https://www.vasp.at/)                          | **vasp**            | r      | 'vasp-out'      | 'OUTCAR'                                                 |
| [WIEN2k](https://susi.theochem.tuwien.ac.at/)         | **lapw1**           | w      | 'lapw-kpt'      | 'case.klist_band'                                        |
| [WIEN2k](https://susi.theochem.tuwien.ac.at/)         | **lapw1**           | r      | 'lapw-out'      | 'case.output1'                                           |
| [BoltzTraP](https://goo.gl/atsFQ8)                    | **BoltzTraP**       | w      | 'boltztrap-in'  | 'case.def', 'case.intrans', 'case.struct', 'case.energy' |
| [BoltzTraP](https://goo.gl/atsFQ8)                    | **BoltzTraP**       | r      | 'boltztrap-out' | 'case.intrans', 'case.trace'                             |
| [BoltzTraP](https://goo.gl/atsFQ8)                    | **BoltzTraP**       | r      | 'boltztrap-dos' | 'case.intrans', 'case.transdos'                          |

#### Some codes use separate files for different spin components

| Code                                                  | Executable      | Access | fileformat      | filenames                                                  |
|-------------------------------------------------------|-----------------|--------|-----------------|------------------------------------------------------------|
| [Quantum Espresso](https://www.quantum-espresso.org/) | **bands.x**     | r      | 'bands-out'     | 'bands.outup', 'bands.outdn'                               |
| [Wannier90](https://www.wannier.org/)                 | **wannier90.x** | r      | 'wannier-out'   | 'seedname_band.datup', 'seedname_band.datdn'               |
| [WIEN2k](https://susi.theochem.tuwien.ac.at/)         | **lapw1**       | r      | 'lapw-out'      | 'case.output1up', 'case.output1dn'                         |
| [BoltzTraP](https://goo.gl/atsFQ8)                    | **BoltzTraP**   | w      | 'boltztrap-in'  | 'case.def', 'case.intrans', 'case.struct', 'case.energyso' |

#### Notes on different codes

[Quantum Espresso](https://www.quantum-espresso.org/)

* **pw.x** must be run with 'verbosity = "high"' to force writing symmetry operations and k-points to file 'prefix.out'.

[Wannier90](https://www.wannier.org/)

* Attribute kpoint of class Kpoint is not available from the input or output files of 'wannier90.x'. It can be constructed from attributes kpath and kindex read from file 'seedname.win' using methods calc_kindex and calc_kpoint.

[VASP](https://www.vasp.at/)

* Symmetry operations are not available from the input or output files of **vasp**. They can be constructed from file 'POSCAR' using [ASE](https://wiki.fysik.dtu.dk/ase/) and [SPGLIB](https://atztogo.github.io/spglib/python-spglib.html) as done in [BoltzTraP](https://goo.gl/atsFQ8) 'vasp2boltz.py'. This is currently not implemented in BRAVE.

[WIEN2k](https://susi.theochem.tuwien.ac.at/)

* File 'case.struct' contains symmetry operations in cartesian coordinates. They can be read by BRAVE and converted to crystal coordinates. This is currently not implemented in BRAVE.
* If **lapw1** is run in parallel file 'case.output1' can be gathered by running 'spaghetti' or manually
```
cat case.output1_? > case.output1
cat case.output1_?? >> case.output1
```
