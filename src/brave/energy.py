"""This module defines class Energy."""

import numpy as np

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
        if not isinstance(value, np.ndarray):
            raise TypeError('energy {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 3:
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
        """Sets the new value of eunit and converts energy, efermi and vref.

    Args:
        eunit (str): New value of eunit. Possible values are 'ev', 'rydberg',
            'hartree', 'thz' and 'cm-1'.
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

    def calc_efermi(self, soc = None):
        """Calculates efermi from nelec for insulators.

    Args:
        soc (bool): Set to True if the calculation includes the spin-orbit
            coupling.

    Returns:
        evbm (float): VBM energy in units of eunit.
        ecbm (float): CBM energy in units of eunit.
        kvbm (ndarray): VBM k-point in units of kunit.
        kcbm (ndarray): CBM k-point in units of kunit.

    Places efermi in the middle of the band gap between the VBM (valence band
    maximum) and the CBM (conduction band minimum).
        """
        if soc is None:
            soc = False

        if soc:
            spin_degeneracy = 1
        else:
            spin_degeneracy = 2

        nval = self.nelec / spin_degeneracy

        if nval.is_integer():
            nval = int(nval)
        else:
            raise ValueError(nval)

        nkpoint = self.nkpoint
        nband = self.nband
        nspin = self.nspin
        energy = self.energy
        kpoint = self.kpoint
        ivbm = np.unravel_index(energy[:, :nval, :].argmax(), (
                nkpoint, nval, nspin))
        icbm = np.unravel_index(energy[:, nval:, :].argmin(), (
                nkpoint, nband - nval, nspin))
        evbm = energy[ivbm[0], ivbm[1], ivbm[2]]
        ecbm = energy[icbm[0], icbm[1] + nval, icbm[2]]
        kvbm = kpoint[ivbm[0], :]
        kcbm = kpoint[icbm[0], :]

        self.efermi = (evbm + ecbm) / 2.0
        return evbm, ecbm, kvbm, kcbm

    def sort_energy(self):
        """Sorts energy in ascending order by band.

    This fixes discontinuities in the band structure.
        """
        _energy = np.copy(self.energy)

        slice = np.zeros(self.nband, float)
        for ikpoint in range(1, self.nkpoint):
            for ispin in range(self.nspin):
                slice = _energy[ikpoint, :, ispin]
                dummy = np.sort(slice)
                _energy[ikpoint, :, ispin] = dummy

        self.energy = _energy

    def read(self, fileformat, filenames, etype = None, lapwkunit = None):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.
        etype (str): Energy type. Possible values are below.
        lapwkunit (str): WIEN2k workaround. WIEN2k requires k-points in crystal
            coordinates with respect to conventional reciprocal lattice vectors
            (see xcrysden/tests/supportInfo.kpath). Possible values are below.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-out'         ['prefix.out']
    'bands-out'      ['prefix.out', 'bands.out'] or ['prefix.out',
                         'bands.outup', 'bands.outdn']
    'matdyn-out'     ['prefix.out', 'matdyn.modes']
    'inteqp-out'     ['prefix.out', 'bandstructure.dat']
    'sigma-out'      ['prefix.out', 'sigma_hp.log']
    'wannier-out'    ['seedname.win', 'seedname_band.dat'] or ['seedname.win',
                         'seedname_band.datup', 'seedname_band.datdn']
    'vasp-out'       ['OUTCAR']
    'lapw-out'       ['case.output1'] or ['case.output1up', 'case.output1dn']

    fileformat       etype
    ----------       -----
    'inteqp-out'     'emf' or 'eqp'
    'sigma-out'      'edft', 'ecor', 'eqp0', 'eqp1', 'eqp0p' or 'eqp1p'

    fileformat       lapwkunit
    ----------       ---------
    'lapw-out'       'cartesian' (for fcc or bcc) or 'crystal' (for hcp)
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

    def write(self, fileformat, filenames, lapwkunit = None, boltzparam = None):
        """Writes properties to files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.
        lapwkunit (str): WIEN2k workaround. WIEN2k requires k-points in crystal
            coordinates with respect to conventional reciprocal lattice vectors
            (see xcrysden/tests/supportInfo.kpath). Possible values are below.
        boltzparam (list): BoltzTraP parameters. Possible values are below.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-in'          ['prefix.in']
    'wannier-in'     ['seedname.win']
    'vasp-kpt'       ['KPOINTS']
    'lapw-kpt'       ['case.klist_band']
    'boltztrap-in'   ['case.def', 'case.intrans', 'case.struct', 'case.energy']
                         (for spin-unpolarized case) or ['case.def',
                         'case.intrans', 'case.struct', 'case.energyso'] (for
                         spin-polarized case)

    fileformat       lapwkunit
    ----------       ---------
    'lapw-kpt'       'cartesian' (for fcc or bcc) or 'crystal' (for hcp)

    fileformat       boltzparam
    ----------       ----------
    'boltztrap-in'   [0.0005, 0.6, 5, 0.3, 1200.0, 10.0, -1.0, 'TETRA', 0]
        """
        if lapwkunit is None:
            lapwkunit = 'cartesian'
        if boltzparam is None:
            boltzparam = [0.0005, 0.6, 5, 0.3, 1200.0, 10.0, -1.0, 'TETRA', 0]

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
            self._write_file_boltztrap_in(3, filenames, boltzparam)
        else:
            raise ValueError(fileformat)

    def __init__(
            self, eunit = None, energy = None, efermi = None, vref = None,
            **kwargs):
        super().__init__(**kwargs)

        if eunit is not None:
            self.eunit = eunit
        if energy is not None:
            self.energy = energy
        if efermi is not None:
            self.efermi = efermi
        if vref is not None:
            self.vref = vref

