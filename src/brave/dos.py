"""This module defines class DOS."""

import numpy

import brave.common as common
from brave.cell import Cell

class DOS(Cell):
    """Class for representing electron and phonon density of states
    (DOS) as functions of electron and phonon energy.
    """

    @property
    def dunit(self):
        """Unit for dos, list of two strs, ['uc'|'bohr3'|'angstrom3'
    |'nm3', 'ev'|'rydberg'|'hartree'|'thz'|'cm-1'].
        """
        return self._dunit

    @dunit.setter
    def dunit(self, dunit):
        self._dunit = [dunit[0].lower(), dunit[1].lower()]

        if self._dunit[0] not in [
                'uc', 'bohr3', 'angstrom3', 'nm3'] or self._dunit[1] not in [
                'ev', 'rydberg', 'hartree', 'thz', 'cm-1']:
            raise ValueError(dunit)

    @dunit.deleter
    def dunit(self):
        del self._dunit

    @property
    def ndos(self):
        """Number of points for the energy grid, int."""
        return self._dos.shape[1]

    @property
    def dos(self):
        """Density of states, array of 2 by ndos floats,
    in units of dunit[1] and el/(dunit[0] dunit[1])
    (electrons per dunit[0] per dunit[1]) or
    ph/(dunit[0] dunit[1]) (phonons per dunit[0]
    per dunit[1]).
        """
        return self._dos

    @dos.setter
    def dos(self, dos):
        self._dos = numpy.array(dos, float)

        if len(self._dos.shape) != 2 or self._dos.shape[0] != 2:
            raise ValueError(dos)

    @dos.deleter
    def dos(self):
        del self._dos

    def set_dunit(self, dunit):
        """Method for setting the new value of dunit = ['uc'|'bohr3'
    |'angstrom3'|'nm3', 'ev'|'rydberg'|'hartree'|'thz'|'cm-1']
    and recalculating dos.
        """
        olddunit = self.dunit
        self.dunit = dunit

        if self.dunit[0] != olddunit[0]:
            oldaunit = self.aunit
            self.set_aunit('angstrom')
            _dscale = {
                    'uc': 1.0 / (self.avol * self.alat ** 3),
                    'bohr3': common._ascale['bohr'] ** 3,
                    'angstrom3': common._ascale['angstrom'] ** 3,
                    'nm3': common._ascale['nm'] ** 3}
            self.set_aunit(oldaunit)

            if hasattr(self, 'dos'):
                dummy = self.dos
                dummy[1] *= _dscale[olddunit[0]] / _dscale[self.dunit[0]]
                self.dos = dummy

        if self.dunit[1] != olddunit[1]:
            if hasattr(self, 'dos'):
                dummy = self.dos
                dummy[0] *= common._escale[self.dunit[1]] / common._escale[
                        olddunit[1]]
                dummy[1] /= common._escale[self.dunit[1]] / common._escale[
                        olddunit[1]]
                self.dos = dummy

    def read(self, fileformat, filenames):
        """Method for reading properties from file.

    fileformat         filenames
    ----------         ---------
    'boltztrap-dos'    ['case.intrans', 'case.transdos']
    'matdyn-dos'       ['prefix.vdos']

    Inherits fileformat and filenames from class Cell.
        """

        if fileformat.lower() == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames)
        elif fileformat.lower() == 'matdyn-dos':
            self._read_dos_matdyn_dos(filenames)
        else:
            super().read(fileformat, filenames)

    def __init__(self, dunit=None, dos=None, **kwargs):
        super().__init__(**kwargs)

        if dunit != None:
            self.dunit = dunit
        if dos != None:
            self.dos = dos

    def _read_dos_boltztrap_dos(self, filenames):
        contents = common._read_file(filenames)

        nspin = 2
        efermi = float(contents[0][2].split()[0])

        nn = len(contents[1]) - 1
        dos = numpy.empty((2, nn), float)
        for ii in range(nn):
            tt = contents[1][ii + 1].split()
            dos[0, ii] = float(tt[0]) - efermi
            dos[1, ii] = float(tt[1]) * nspin

        self.dunit, self.dos = ['uc', 'rydberg'], dos

    def _read_dos_matdyn_dos(self, filenames):
        contents = common._read_file(filenames)

        nn = len(contents[0])
        dos = numpy.empty((2, nn), float)
        for ii in range(nn):
            tt = contents[0][ii].split()
            dos[0, ii] = float(tt[0])
            dos[1, ii] = float(tt[1])

        self.dunit, self.dos = ['uc', 'cm-1'], dos

