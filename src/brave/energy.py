"""This module defines class Energy."""

import numpy

import brave.common as common
from brave.kpoint import Kpoint

class Energy(Kpoint):
    """Class for representing single-particle energies of electrons
    and phonons. Used for converting electron and phonon dispersion
    relations to formats suitable for various post-processing codes.
    """

    @property
    def eunit(self):
        """Unit for energy, str, 'ev'|'rydberg'|'hartree'|'thz'
    |'cm-1'.
        """
        return self._eunit

    @eunit.setter
    def eunit(self, eunit):
        self._eunit = eunit.lower()

        if self._eunit not in ['ev', 'rydberg', 'hartree', 'thz', 'cm-1']:
            raise ValueError(eunit)

    @eunit.deleter
    def eunit(self):
        del self._eunit

    @property
    def nband(self):
        """Number of electron bands or phonon modes, int."""
        return self._energy.shape[1]

    @property
    def nspin(self):
        """Number of electron spins or 1 for phonons, int."""
        return self._energy.shape[2]

    @property
    def energy(self):
        """Energy values, array of nkpoint by nband by nspin floats,
    in units of eunit.
        """
        return self._energy

    @energy.setter
    def energy(self, energy):
        self._energy = numpy.array(energy, float)

        _nkpoint = -1
        if hasattr(self, 'nkpoint'):
            _nkpoint = self.nkpoint

        if len(self._energy.shape) != 3 or self._energy.shape[0] != _nkpoint:
            raise ValueError(energy)

    @energy.deleter
    def energy(self):
        del self._energy

    @property
    def efermi(self):
        """Fermi level, float, in units of eunit.
        """
        return self._efermi

    @efermi.setter
    def efermi(self, efermi):
        self._efermi = float(efermi)

    @efermi.deleter
    def efermi(self):
        del self._efermi

    @property
    def vref(self):
        """Reference potential, float, in units of eunit.
        """
        return self._vref

    @vref.setter
    def vref(self, vref):
        self._vref = float(vref)

    @vref.deleter
    def vref(self):
        del self._vref

    def set_eunit(self, eunit):
        """Method for setting the new value of eunit =
    'ev'|'rydberg'|'hartree'|'thz'|'cm-1' and
    recalculating energy, efermi and vref.
        """
        oldeunit = self.eunit
        self.eunit = eunit

        if self.eunit != oldeunit:
            if hasattr(self, 'energy'):
                self.energy *= common._escale[self.eunit] / common._escale[
                        oldeunit]

            if hasattr(self, 'efermi'):
                self.efermi *= common._escale[self.eunit] / common._escale[
                        oldeunit]

            if hasattr(self, 'vref'):
                self.vref *= common._escale[self.eunit] / common._escale[
                        oldeunit]

    def calc_efermi(self):
        """Method for calculating efermi for insulators given the
    number of electrons, nelec. Places efermi in the
    middle of the band gap. Returns energies evbm and
    ecbm in units of eunit and k-points kvbm and kcbm
    in units of kunit at the vbm (valence band maximum)
    and cbm (conduction band minimum).
        """
        if round(self.nelec) != self.nelec or int(round(self.nelec)) % 2 != 0:
            raise ValueError(self.nelec)
        nval = int(round(self.nelec)) // 2

        nkpoint = self.nkpoint
        nband = self.nband
        nspin = self.nspin
        energy = self.energy
        kpoint = self.kpoint
        _ivbm = numpy.unravel_index(numpy.argmax(energy[:, :nval, :]), (
                nkpoint, nval, nspin))
        _icbm = numpy.unravel_index(numpy.argmin(energy[:, nval:, :]), (
                nkpoint, nband - nval, nspin))
        _evbm = energy[_ivbm[0], _ivbm[1], _ivbm[2]]
        _ecbm = energy[_icbm[0], _icbm[1] + nval, _icbm[2]]
        _kvbm = kpoint[_ivbm[0]]
        _kcbm = kpoint[_icbm[0]]

        self.efermi = (_evbm + _ecbm) / 2.0
        return _evbm, _ecbm, _kvbm, _kcbm

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

    def read(self, fileformat, filenames, etype=None, lapwkunit='cartesian'):
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

        if fileformat.lower() == 'internal':
            self._read_file_internal(3, filenames)
        elif fileformat.lower() == 'pw-out':
            self._read_file_pw_out(3, filenames)
        elif fileformat.lower() == 'bands-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_bands_out(3, filenames[1:])
        elif fileformat.lower() == 'matdyn-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_matdyn_out(3, [filenames[1]])
        elif fileformat.lower() == 'inteqp-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_inteqp_out(3, [filenames[1]], etype)
        elif fileformat.lower() == 'sigma-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_sigma_out(3, [filenames[1]], etype)
        elif fileformat.lower() == 'wannier-out':
            self._read_file_wannier_in(1, [filenames[0]])
            self._read_file_wannier_out(3, filenames[1:])
        elif fileformat.lower() == 'vasp-out':
            self._read_file_vasp_out(3, filenames)
        elif fileformat.lower() == 'lapw-out':
            self._read_file_lapw_out(3, filenames, lapwkunit)
        else:
            raise ValueError(fileformat)

    def write(
            self, fileformat, filenames, lapwkunit='cartesian',
            deltae=0.0005, ecut=0.6, lpfac=5, efcut=0.3, tmax=1200.0,
            deltat=10.0, ecut2=-1.0, dosmethod='TETRA', nband_exclude=0):
        """Method for writing properties to file.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-in'          ['prefix.in']
    'wannier-in'     ['seedname.win']
    'vasp-kpt'       ['KPOINTS']
    'lapw-kpt'       ['case.klist_band']
    'boltztrap-in'   ['prefix.def', 'prefix.intrans',
                             'prefix.struct', 'prefix.energy[so]']

    for 'boltztrap-in' use 'prefix.energy' for spin-unpolarized
            case and 'prefix.energyso' for spin-polarized or
            non-collinear case

    lapwkunit = 'cartesian'|'crystal' used for 'lapw-kpt',
            WIEN2k requires k-points in crystal coordinates
            with respect to conventional reciprocal lattice
            vectors, see xcrysden/tests/supportInfo.kpath
            for details, use 'cartesian' for fcc or bcc
            and 'crystal' for hcp
    deltae, ecut, lpfac, efcut, tmax, deltat, ecut2, dosmethod,
            nband_exclude used for 'boltztrap-in'
        """

        if fileformat.lower() == 'internal':
            self._write_file_internal(3, filenames)
        elif fileformat.lower() == 'pw-in':
            self._write_file_pw_in(2, filenames)
        elif fileformat.lower() == 'wannier-in':
            self._write_file_wannier_in(2, filenames)
        elif fileformat.lower() == 'vasp-kpt':
            self._write_file_vasp_kpt(2, filenames)
        elif fileformat.lower() == 'lapw-kpt':
            self._write_file_lapw_kpt(2, filenames, lapwkunit) 
        elif fileformat.lower() == 'boltztrap-in':
            self._write_file_boltztrap_in(
                    3, filenames, deltae, ecut, lpfac, efcut, tmax, deltat,
                    ecut2, dosmethod, nband_exclude)
        else:
            raise ValueError(fileformat)

    def __init__(
            self, eunit=None, energy=None, efermi=None, vref=None, **kwargs):
        super().__init__(**kwargs)

        if eunit != None:
            self.eunit = eunit
        if energy != None:
            self.energy = energy
        if efermi != None:
            self.efermi = efermi
        if vref != None:
            self.vref = vref

