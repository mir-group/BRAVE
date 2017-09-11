"""This module defines class Cell."""

import numpy

import brave.common as common
from brave.file import File

class Cell(File):
    """Class for representing the lattice structure.

    Class Cell defines the direct and reciprocal lattice vectors and the
    symmetry operations. It is used for converting between Cartesian and
    crystal coordinates in class Kpoint and between different units in classes
    DOS and Transport.
    """

    @property
    def prefix(self):
        """A string holding the name of the system, same as prefix in Quantum
    ESPRESSO, seedname in Wannier90, SYSTEM in VASP, or case in WIEN2k.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, value):
        if not isinstance(value, str):
            raise TypeError('prefix {0!r}'.format(value))
        self._prefix = value

    @prefix.deleter
    def prefix(self):
        del self._prefix

    @property
    def aunit(self):
        """A string holding the units of alat. Possible values are 'bohr',
    'angstrom' and 'nm'.
        """
        return self._aunit

    @aunit.setter
    def aunit(self, value):
        if not isinstance(value, str):
            raise TypeError('aunit {0!r}'.format(value))
        if value not in common._ascale.keys():
            raise ValueError('aunit {0!r}'.format(value))
        self._aunit = value

    @aunit.deleter
    def aunit(self):
        del self._aunit

    @property
    def alat(self):
        """A float holding the direct lattice constant, in units of aunit."""
        return self._alat

    @alat.setter
    def alat(self, value):
        if not isinstance(value, float):
            raise TypeError('alat {0!r}'.format(value))
        self._alat = value

    @alat.deleter
    def alat(self):
        del self._alat

    @property
    def avec(self):
        """A 3 by 3 ndarray of floats holding the direct lattice vectors, in
    units of alat.
        """
        return self._avec

    @avec.setter
    def avec(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('avec {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or value.shape != (3, 3):
            raise ValueError('avec {0!r}'.format(value))
        self._avec = value

    @avec.deleter
    def avec(self):
        del self._avec

    @property
    def bvec(self):
        """A 3 by 3 ndarray of floats holding the reciprocal lattice vectors,
    in units of 2*pi/alat.
        """
        return self._bvec

    @bvec.setter
    def bvec(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('bvec {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or value.shape != (3, 3):
            raise ValueError('bvec {0!r}'.format(value))
        self._bvec = value

    @bvec.deleter
    def bvec(self):
        del self._bvec

    @property
    def avol(self):
        """A float holding the volume of the direct lattice primitive cell, in
    units of alat^3.
        """
        return self._avol

    @avol.setter
    def avol(self, value):
        if not isinstance(value, float):
            raise TypeError('avol {0!r}'.format(value))
        self._avol = value

    @avol.deleter
    def avol(self):
        del self._avol

    @property
    def bvol(self):
        """A float holding the volume of the reciprocal lattice primitive cell,
    in units of (2*pi/alat)^3.
        """
        return self._bvol

    @bvol.setter
    def bvol(self, value):
        if not isinstance(value, float):
            raise TypeError('bvol {0!r}'.format(value))
        self._bvol = value

    @bvol.deleter
    def bvol(self):
        del self._bvol

    @property
    def natom(self):
        """An integer holding the number of atoms per primitive cell."""
        return self._natom

    @natom.setter
    def natom(self, value):
        if not isinstance(value, int):
            raise TypeError('natom {0!r}'.format(value))
        self._natom = value

    @natom.deleter
    def natom(self):
        del self._natom

    @property
    def nelec(self):
        """A float holding the number of electrons per primitive cell."""
        return self._nelec

    @nelec.setter
    def nelec(self, value):
        if not isinstance(value, float):
            raise TypeError('nelec {0!r}'.format(value))
        self._nelec = value

    @nelec.deleter
    def nelec(self):
        del self._nelec

    @property
    def nsym(self):
        """An integer holding the number of symmetry operations."""
        return self._rot.shape[0]

    @property
    def rot(self):
        """An nsym by 3 by 3 ndarray of integers holding the rotation matrices
    operating on crystal coordinates.
        """
        return self._rot

    @rot.setter
    def rot(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('rot {0!r}'.format(value))
        if value.dtype != numpy.dtype('int') or len(
                value.shape) != 3 or value.shape[1:3] != (3, 3):
            raise ValueError('rot {0!r}'.format(value))
        self._rot = value

    @rot.deleter
    def rot(self):
        del self._rot

    def set_aunit(self, aunit):
        """Method for setting the new value of aunit =
    'bohr'|'angstrom'|'nm' and recalculating alat.
        """
        if aunit != self.aunit:
            if hasattr(self, 'alat'):
                self.alat *= common._ascale[aunit] / common._ascale[self.aunit]
            self.aunit = aunit

    def set_alat(self, alat):
        """Method for setting the new value of alat = float
    and recalculating avec, bvec, avol, and bvol.
        """
        if abs(alat - self.alat) > common.EPS12:
            if hasattr(self, 'avec'):
                self.avec *= self.alat / alat
            if hasattr(self, 'bvec'):
                self.bvec *= alat / self.alat
            if hasattr(self, 'avol'):
                self.avol *= (self.alat / alat) ** 3
            if hasattr(self, 'bvol'):
                self.bvol *= (alat / self.alat) ** 3
            self.alat = alat

    def calc_avec(self):
        """Method for calculating avec given bvec."""
        self.avec = self._calc_vec(self.bvec)

    def calc_bvec(self):
        """Method for calculating bvec given avec."""
        self.bvec = self._calc_vec(self.avec)

    def calc_avol(self):
        """Method for calculating avol given avec."""
        self.avol = abs(self._calc_vol(self.avec))

    def calc_bvol(self):
        """Method for calculating bvol given bvec."""
        self.bvol = abs(self._calc_vol(self.bvec))

    def read(self, fileformat, filenames):
        """Method for reading properties from file.

    fileformat      filenames
    ----------      ---------
    'internal'      ['prefix.brave']
    'pw-out'        ['prefix.out']
    'wannier-in'    ['seedname.win']
    'vasp-out'      ['OUTCAR']
    'lapw-out'      ['case.output1']
                 or ['case.output1up', 'case.output1dn']
        """

        if fileformat == 'internal':
            self._read_file_internal(1, filenames)
        elif fileformat == 'pw-out':
            self._read_file_pw_out(1, filenames)
        elif fileformat == 'wannier-in':
            self._read_file_wannier_in(1, filenames)
        elif fileformat == 'vasp-out':
            self._read_file_vasp_out(1, filenames)
        elif fileformat == 'lapw-out':
            self._read_file_lapw_out(1, filenames, None)
        else:
            raise ValueError(fileformat)

    def write(self, fileformat, filenames):
        """Method for writing properties to file.

    fileformat      filenames
    ----------      ---------
    'internal'      ['prefix.brave']
    'pw-in'         ['prefix.in']
    'wannier-in'    ['seedname.win']
        """

        if fileformat == 'internal':
            self._write_file_internal(1, filenames)
        elif fileformat == 'pw-in':
            self._write_file_pw_in(1, filenames)
        elif fileformat == 'wannier-in':
            self._write_file_wannier_in(1, filenames)
        else:
            raise ValueError(fileformat)

    def _calc_vec(self, xvec):
        yvec = numpy.empty((3, 3), float)
        yvec[0] = numpy.cross(xvec[1], xvec[2])
        yvec[1] = numpy.cross(xvec[2], xvec[0])
        yvec[2] = numpy.cross(xvec[0], xvec[1])
        yvec /= self._calc_vol(xvec)
        return yvec

    def _calc_vol(self, vec):
        vol = numpy.dot(vec[0], numpy.cross(vec[1], vec[2]))
        return vol

    def __init__(
            self, prefix=None, aunit=None, alat=None, avec=None, bvec=None,
            avol=None, bvol=None, natom=None, nelec=None, rot=None):

        if prefix is not None:
            self.prefix = prefix
        if aunit is not None:
            self.aunit = aunit
        if alat is not None:
            self.alat = alat
        if avec is not None:
            self.avec = avec
        if bvec is not None:
            self.bvec = bvec
        if avol is not None:
            self.avol = avol
        if bvol is not None:
            self.bvol = bvol
        if natom is not None:
            self.natom = natom
        if nelec is not None:
            self.nelec = nelec
        if rot is None:
            self.rot = numpy.array([numpy.identity(3, int)])
        else:
            self.rot = rot

