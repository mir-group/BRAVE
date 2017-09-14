"""This module defines class Energy."""

import numpy

import brave.common as common
from brave.kpoint import Kpoint

class Energy(Kpoint):
    """Class for representing the energy bands.

    Class Energy defines the electron or phonon energy bands. It is used for
    plotting the energy band diagrams in class Diagram and for generating the
    input files for calculating the electronic transport coefficients.
    """

    @property
    def eunit(self):
        """A string holding the units of energy. Possible values are 'ev',
        'rydberg', 'hartree', 'thz' and 'cm-1'.
        """
        return self._eunit

    @eunit.setter
    def eunit(self, value):
        if not isinstance(value, str):
            raise TypeError('eunit {0!r}'.format(value))
        if value not in common._escale.keys():
            raise ValueError('eunit {0!r}'.format(value))
        self._eunit = value

    @eunit.deleter
    def eunit(self):
        del self._eunit

    @property
    def nkpoint(self):
        """An integer holding the number of k-points."""
        _list = []
        if hasattr(self, 'energy'):
            _list.append(self._energy.shape[0])
        if hasattr(self, 'kpoint') or hasattr(self, 'kline') or hasattr(
                self, 'kweight'):
            _list.append(super().nkpoint)

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nkpoint')

    @property
    def nband(self):
        """An integer holding the number of electron energy bands or phonon
    modes.
        """
        return self._energy.shape[1]

    @property
    def nspin(self):
        """An integer holding the number of spin components for electrons or 1
    for phonons.
        """
        return self._energy.shape[2]

    @property
    def energy(self):
        """A nkpoint by nband by nspin ndarray of floats holding the electron
    or phonon energies, in units of eunit.
        """
        return self._energy

    @energy.setter
    def energy(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('energy {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 3:
            raise ValueError('energy {0!r}'.format(value))
        self._energy = value

    @energy.deleter
    def energy(self):
        del self._energy

    @property
    def efermi(self):
        """A float holding the chemical potential of electrons or the Fermi
    level, in units of eunit.
        """
        return self._efermi

    @efermi.setter
    def efermi(self, value):
        if not isinstance(value, float):
            raise TypeError('efermi {0!r}'.format(value))
        self._efermi = value

    @efermi.deleter
    def efermi(self):
        del self._efermi

    @property
    def vref(self):
        """A float holding the reference potential, in units of eunit."""
        return self._vref

    @vref.setter
    def vref(self, value):
        if not isinstance(value, float):
            raise TypeError('vref {0!r}'.format(value))
        self._vref = value

    @vref.deleter
    def vref(self):
        del self._vref

    def set_eunit(self, eunit):
        """Method for setting the new value of eunit =
    'ev'|'rydberg'|'hartree'|'thz'|'cm-1' and
    recalculating energy, efermi and vref.
        """
        if eunit != self.eunit:
            if hasattr(self, 'energy'):
                self.energy *= common._escale[eunit] / common._escale[
                        self.eunit]
            if hasattr(self, 'efermi'):
                self.efermi *= common._escale[eunit] / common._escale[
                        self.eunit]
            if hasattr(self, 'vref'):
                self.vref *= common._escale[eunit] / common._escale[self.eunit]
            self.eunit = eunit

    def calc_efermi(self, soc=None):
        """Method for calculating efermi for insulators given the
    number of electrons, nelec. Places efermi in the
    middle of the band gap. Returns energies evbm and
    ecbm in units of eunit and k-points kvbm and kcbm
    in units of kunit at the vbm (valence band maximum)
    and cbm (conduction band minimum). Set soc to True
    if the calculation includes the spin-orbit coupling.
        """
        if soc is None:
            soc = False

        nval = self.nelec
        if not soc:
            nval /= 2
        if nval.is_integer():
            nval = int(nval)
        else:
            raise ValueError(nval)

        nkpoint = self.nkpoint
        nband = self.nband
        nspin = self.nspin
        energy = self.energy
        kpoint = self.kpoint
        ivbm = numpy.unravel_index(energy[:, :nval, :].argmax(), (
                nkpoint, nval, nspin))
        icbm = numpy.unravel_index(energy[:, nval:, :].argmin(), (
                nkpoint, nband - nval, nspin))
        evbm = energy[ivbm[0], ivbm[1], ivbm[2]]
        ecbm = energy[icbm[0], icbm[1] + nval, icbm[2]]
        kvbm = kpoint[ivbm[0], :]
        kcbm = kpoint[icbm[0], :]

        self.efermi = (evbm + ecbm) / 2.0
        return evbm, ecbm, kvbm, kcbm

    def sort_energy(self):
        """Method for sorting energy in ascending order by band
    to fix discontinuities in the band structure.
        """
        _energy = numpy.copy(self.energy)

        slice = numpy.zeros(self.nband, float)
        for ikpoint in range(1, self.nkpoint):
            for ispin in range(self.nspin):
                slice = _energy[ikpoint, :, ispin]
                dummy = numpy.sort(slice)
                _energy[ikpoint, :, ispin] = dummy

        self.energy = _energy

    def read(self, fileformat, filenames, etype=None, lapwkunit=None):
        """Method for reading properties from file.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-out'         ['prefix.out']
    'bands-out'      ['prefix.out', 'bands.out']
                  or ['prefix.out', 'bands.outup', 'bands.outdn']
    'matdyn-out'     ['prefix.out', 'matdyn.modes']
    'inteqp-out'     ['prefix.out', 'bandstructure.dat']
    'sigma-out'      ['prefix.out', 'sigma_hp.log']
    'wannier-out'    ['seedname.win', 'seedname_band.dat']
                  or ['seedname.win', 'seedname_band.datup',
                             'seedname_band.datdn']
    'vasp-out'       ['OUTCAR']
    'lapw-out'       ['case.output1']
                  or ['case.output1up', 'case.output1dn']

    etype = 'emf'|'eqp' used for 'inteqp-out'
    etype = 'edft'|'ecor'|'eqp0'|'eqp1'|'eqp0p'|'eqp1p' used
            for 'sigma-out'
    lapwkunit = 'cartesian'|'crystal' used for 'lapw-kpt',
            WIEN2k requires k-points in crystal coordinates
            with respect to conventional reciprocal lattice
            vectors, see xcrysden/tests/supportInfo.kpath
            for details, use 'cartesian' for fcc or bcc
            and 'crystal' for hcp
        """
        if lapwkunit is None:
            lapwkunit = 'cartesian'

        if fileformat == 'internal':
            self._read_file_internal(3, filenames)
        elif fileformat == 'pw-out':
            self._read_file_pw_out(3, filenames)
        elif fileformat == 'bands-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_bands_out(3, filenames[1:])
        elif fileformat == 'matdyn-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_matdyn_out(3, [filenames[1]])
        elif fileformat == 'inteqp-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_inteqp_out(3, [filenames[1]], etype)
        elif fileformat == 'sigma-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_sigma_out(3, [filenames[1]], etype)
        elif fileformat == 'wannier-out':
            self._read_file_wannier_in(1, [filenames[0]])
            self._read_file_wannier_out(3, filenames[1:])
        elif fileformat == 'vasp-out':
            self._read_file_vasp_out(3, filenames)
        elif fileformat == 'lapw-out':
            self._read_file_lapw_out(3, filenames, lapwkunit)
        else:
            raise ValueError(fileformat)

    def write(
            self, fileformat, filenames, lapwkunit=None, deltae=None,
            ecut=None, lpfac=None, efcut=None, tmax=None, deltat=None,
            ecut2=None, dosmethod=None, nband_exclude=None):
        """Method for writing properties to file.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-in'          ['prefix.in']
    'wannier-in'     ['seedname.win']
    'vasp-kpt'       ['KPOINTS']
    'lapw-kpt'       ['case.klist_band']
    'boltztrap-in'   ['case.def', 'case.intrans', 'case.struct',
                             'case.energy[so]']

    for 'boltztrap-in' use 'case.energy' for spin-unpolarized
            case and 'case.energyso' for spin-polarized case

    lapwkunit = 'cartesian'|'crystal' used for 'lapw-kpt',
            WIEN2k requires k-points in crystal coordinates
            with respect to conventional reciprocal lattice
            vectors, see xcrysden/tests/supportInfo.kpath
            for details, use 'cartesian' for fcc or bcc
            and 'crystal' for hcp
    deltae, ecut, lpfac, efcut, tmax, deltat, ecut2, dosmethod,
            nband_exclude used for 'boltztrap-in'
        """
        if lapwkunit is None:
            lapwkunit = 'cartesian'
        if deltae is None:
            deltae = 0.0005
        if ecut is None:
            ecut = 0.6
        if lpfac is None:
            lpfac = 5
        if efcut is None:
            efcut = 0.3
        if tmax is None:
            tmax = 1200.0
        if deltat is None:
            deltat = 10.0
        if ecut2 is None:
            ecut2 = -1.0
        if dosmethod is None:
            dosmethod = 'TETRA'
        if nband_exclude is None:
            nband_exclude = 0

        if fileformat == 'internal':
            self._write_file_internal(3, filenames)
        elif fileformat == 'pw-in':
            self._write_file_pw_in(2, filenames)
        elif fileformat == 'wannier-in':
            self._write_file_wannier_in(2, filenames)
        elif fileformat == 'vasp-kpt':
            self._write_file_vasp_kpt(2, filenames)
        elif fileformat == 'lapw-kpt':
            self._write_file_lapw_kpt(2, filenames, lapwkunit) 
        elif fileformat == 'boltztrap-in':
            self._write_file_boltztrap_in(
                    3, filenames, deltae, ecut, lpfac, efcut, tmax, deltat,
                    ecut2, dosmethod, nband_exclude)
        else:
            raise ValueError(fileformat)

    def __init__(
            self, eunit=None, energy=None, efermi=None, vref=None, **kwargs):
        super().__init__(**kwargs)

        if eunit is not None:
            self.eunit = eunit
        if energy is not None:
            self.energy = energy
        if efermi is not None:
            self.efermi = efermi
        if vref is not None:
            self.vref = vref

