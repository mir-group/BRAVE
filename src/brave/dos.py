"""This module defines class DOS."""

import linecache
import numpy as np

import brave.common as common
from brave.cell import Cell

class DOS(Cell):
    """Class for representing the density of states.

    Class DOS defines the electron or phonon density of states as a function of
    the electron or phonon energy.
    """

    @property
    def dunit(self):
        """A list of 2 strings holding the units of dos. Possible values of the
    first string are 'uc', 'bohr3', 'angstrom3' and 'nm3'. Possible values of
    the second string are 'ev', 'rydberg', 'hartree', 'thz' and 'cm-1'.
        """
        return self._dunit

    @dunit.setter
    def dunit(self, value):
        if not isinstance(value, list):
            raise TypeError('dunit {0!r}'.format(value))
        if len(value) != 2 or value[0] not in common._a3scale.keys() or value[
                1] not in common._escale.keys():
            raise ValueError('dunit {0!r}'.format(value))
        self._dunit = value

    @dunit.deleter
    def dunit(self):
        del self._dunit

    @property
    def ndos(self):
        """An integer holding the number of points in the energy grid."""
        return self._dos.shape[1]

    @property
    def dos(self):
        """A 2 by ndos ndarray of floats holding the density of states per
    spin channel. The first column is in units of dunit[1], the second column
    is in units of 1/(dunit[0] dunit[1]).
        """
        return self._dos

    @dos.setter
    def dos(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('dos {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(
                value.shape) != 2 or value.shape[0] != 2:
            raise ValueError('dos {0!r}'.format(value))
        self._dos = value

    @dos.deleter
    def dos(self):
        del self._dos

    def set_dunit(self, dunit):
        """Sets the new value of dunit and recalculates dos.

    Args:
        dunit (list): New value of dunit. Possible values of dunit[0] are 'uc',
            'bohr3', 'angstrom3' and 'nm3'. Possible values of dunit[1] are
            'ev', 'rydberg', 'hartree', 'thz' and 'cm-1'.
        """
        if dunit[0] != self.dunit[0]:
            oldaunit = self.aunit
            self.set_aunit('angstrom')
            common._a3scale['uc'] = 1.0 / (self.avol * self.alat ** 3)
            self.set_aunit(oldaunit)

            if hasattr(self, 'dos'):
                dummy = self.dos
                dummy[1, :] *= common._a3scale[self.dunit[0]] / common._a3scale[
                    dunit[0]]
                self.dos = dummy

        if dunit[1] != self.dunit[1]:
            if hasattr(self, 'dos'):
                dummy = self.dos
                dummy[0, :] *= common._escale[dunit[1]] / common._escale[
                    self.dunit[1]]
                dummy[1, :] /= common._escale[dunit[1]] / common._escale[
                    self.dunit[1]]
                self.dos = dummy

        if dunit != self.dunit:
            self.dunit = dunit

    def read(self, fileformat, filenames, **kwargs):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.

    fileformat         filenames
    ----------         ---------
    'boltztrap-dos'    ['case.intrans', 'case.transdos']
    'matdyn-dos'       ['prefix.vdos']

    Inherits fileformat and filenames from class Cell.
        """
        if fileformat == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames)
        elif fileformat == 'matdyn-dos':
            self._read_dos_matdyn_dos(filenames)
        else:
            super().read(fileformat, filenames, **kwargs)

    def _read_dos_boltztrap_dos(self, filenames):

        efermi = float(linecache.getline(filenames[0], 3).split()[0])
        dos = np.loadtxt(filenames[1], dtype = float, skiprows = 1, usecols = (
                0, 1), unpack = True)
        dos[0, :] -= efermi

        self.dunit, self.dos = ['uc', 'rydberg'], dos

    def _read_dos_matdyn_dos(self, filenames):

        dos = np.loadtxt(filenames[0], dtype = float, skiprows = 1, usecols = (
                0, 1), unpack = True)

        self.dunit, self.dos = ['uc', 'cm-1'], dos

    def __init__(self, dunit=None, dos=None, **kwargs):
        super().__init__(**kwargs)

        if dunit is not None:
            self.dunit = dunit
        if dos is not None:
            self.dos = dos

