"""This module defines class Cell."""

import numpy

import brave.common as common
from brave.file import File

class Cell(File):
    """Class for representing direct and reciprocal lattice vectors
    and symmetry operations. Used for converting between Cartesian
    and crystal coordinates in class Kpoint.
    """

    @property
    def prefix(self):
        """Name of the system, str, same as prefix in Quantum ESPRESSO,
    SYSTEM in VASP, or case in WIEN2k.
        """
        return self._prefix

    @prefix.setter
    def prefix(self, prefix):
        self._prefix = str(prefix)

    @prefix.deleter
    def prefix(self):
        del self._prefix

    @property
    def aunit(self):
        """Unit for alat, str, 'bohr'|'angstrom'|'nm'."""
        return self._aunit

    @aunit.setter
    def aunit(self, aunit):
        self._aunit = aunit.lower()

        if self._aunit not in ['bohr', 'angstrom', 'nm']:
            raise ValueError(aunit)

    @aunit.deleter
    def aunit(self):
        del self._aunit

    @property
    def alat(self):
        """Direct lattice constant, float, in units of aunit."""
        return self._alat

    @alat.setter
    def alat(self, alat):
        self._alat = float(alat)

    @alat.deleter
    def alat(self):
        del self._alat

    @property
    def avec(self):
        """Direct lattice vectors, array of 3 by 3 floats,
    in units of alat.
        """
        return self._avec

    @avec.setter
    def avec(self, avec):
        self._avec = numpy.array(avec, float)

        if self._avec.shape != (3, 3):
            raise ValueError(avec)

    @avec.deleter
    def avec(self):
        del self._avec

    @property
    def bvec(self):
        """Reciprocal lattice vectors, array of 3 by 3 floats,
    in units of 2*pi/alat.
        """
        return self._bvec

    @bvec.setter
    def bvec(self, bvec):
        self._bvec = numpy.array(bvec, float)

        if self._bvec.shape != (3, 3):
            raise ValueError(bvec)

    @bvec.deleter
    def bvec(self):
        del self._bvec

    @property
    def avol(self):
        """Volume of direct lattice primitive cell, float,
    in units of alat^3.
        """
        return self._avol

    @avol.setter
    def avol(self, avol):
        self._avol = float(avol)

    @avol.deleter
    def avol(self):
        del self._avol

    @property
    def bvol(self):
        """Volume of reciprocal lattice primitive cell, float,
    in units of (2*pi/alat)^3.
        """
        return self._bvol

    @bvol.setter
    def bvol(self, bvol):
        self._bvol = float(bvol)

    @bvol.deleter
    def bvol(self):
        del self._bvol

    @property
    def natom(self):
        """Number of atoms, int."""
        return self._natom

    @natom.setter
    def natom(self, natom):
        self._natom = int(natom)

    @natom.deleter
    def natom(self):
        del self._natom

    @property
    def nelec(self):
        """Number of electrons, float."""
        return self._nelec

    @nelec.setter
    def nelec(self, nelec):
        self._nelec = float(nelec)

    @nelec.deleter
    def nelec(self):
        del self._nelec

    @property
    def nsym(self):
        """Number of symmetry operations, int."""
        return self._rot.shape[0]

    @property
    def rot(self):
        """Rotational matrices in crystal coordinates, array of
    nsym by 3 by 3 ints.
        """
        return self._rot

    @rot.setter
    def rot(self, rot):
        self._rot = numpy.array(rot, int)

        if len(self._rot.shape) != 3 or self._rot.shape[1] != 3 or (
                self._rot.shape[2] != 3):
            raise ValueError(rot)

    @rot.deleter
    def rot(self):
        del self._rot

    def set_aunit(self, aunit):
        """Method for setting the new value of aunit =
    'bohr'|'angstrom'|'nm' and recalculating alat.
        """
        oldaunit = self.aunit
        self.aunit = aunit

        if self.aunit != oldaunit:
            if hasattr(self, 'alat'):
                self.alat *= common._ascale[
                        self.aunit] / common._ascale[oldaunit]

    def set_alat(self, alat):
        """Method for setting the new value of alat = float
    and recalculating avec, bvec, avol, and bvol.
        """
        oldalat = self.alat
        self.alat = alat

        if abs(self.alat - oldalat) > common.EPS12:
            if hasattr(self, 'avec'):
                self.avec *= oldalat / self.alat

            if hasattr(self, 'bvec'):
                self.bvec *= self.alat / oldalat

            if hasattr(self, 'avol'):
                self.avol *= (oldalat / self.alat) ** 3

            if hasattr(self, 'bvol'):
                self.bvol *= (self.alat / oldalat) ** 3

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

        if fileformat.lower() == 'internal':
            self._read_file_internal(1, filenames)
        elif fileformat.lower() == 'pw-out':
            self._read_file_pw_out(1, filenames)
        elif fileformat.lower() == 'wannier-in':
            self._read_file_wannier_in(1, filenames)
        elif fileformat.lower() == 'vasp-out':
            self._read_file_vasp_out(1, filenames)
        elif fileformat.lower() == 'lapw-out':
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

        if fileformat.lower() == 'internal':
            self._write_file_internal(1, filenames)
        elif fileformat.lower() == 'pw-in':
            self._write_file_pw_in(1, filenames)
        elif fileformat.lower() == 'wannier-in':
            self._write_file_wannier_in(1, filenames)
        else:
            raise ValueError(fileformat)

    def __init__(
            self, prefix=None, aunit=None, alat=None, avec=None, bvec=None,
            avol=None, bvol=None, natom=None, nelec=None, rot=None):

        if prefix != None:
            self.prefix = prefix
        if aunit != None:
            self.aunit = aunit
        if alat != None:
            self.alat = alat
        if avec != None:
            self.avec = avec
        if bvec != None:
            self.bvec = bvec
        if avol != None:
            self.avol = avol
        if bvol != None:
            self.bvol = bvol
        if natom != None:
            self.natom = natom
        if nelec != None:
            self.nelec = nelec
        if rot == None:
            self.rot = numpy.array([numpy.identity(3, float)])
        else:
            self.rot = rot

