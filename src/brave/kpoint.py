"""This module defines class Kpoint."""

import numpy

import brave.common as common
from brave.cell import Cell

class Kpoint(Cell):
    """Class for representing the k-points.

    Class Kpoint defines the grid of k- or q-points in reciprocal space and the
    path along high symmetry lines in reciprocal space. It is used for
    constructing the electron or phonon energy bands in class Energy.
    """

    @property
    def kunit(self):
        """A string holding the units of kpoint and kpath. Possible values are
    'cartesian' and 'crystal', indicating units of 2*pi/alat and bvec,
    respectively.
        """
        return self._kunit

    @kunit.setter
    def kunit(self, kunit):
        self._kunit = kunit

        if self._kunit not in ['cartesian', 'crystal']:
            raise ValueError(kunit)

    @kunit.deleter
    def kunit(self):
        del self._kunit

    @property
    def nkpoint(self):
        """An integer holding the number of k-points."""
        _list = []
        if hasattr(self, 'kpoint'):
            _list.append(self._kpoint.shape[0])
        if hasattr(self, 'kline'):
            _list.append(self._kline.shape[0])
        if hasattr(self, 'kweight'):
            _list.append(self._kweight.shape[0])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nkpoint')

    @property
    def kpoint(self):
        """An nkpoint by 3 ndarray of floats holding the coordinates of
    k-points, in units of kunit.
        """
        return self._kpoint

    @kpoint.setter
    def kpoint(self, kpoint):
        self._kpoint = numpy.array(kpoint, float)

        if len(self._kpoint.shape) != 2 or self._kpoint.shape[1] != 3:
            raise ValueError(kpoint)

    @kpoint.deleter
    def kpoint(self):
        del self._kpoint

    @property
    def kline(self):
        """A length-nkpoint ndarray of floats holding the Cartesian coordinate
    along the path in reciprocal space, in units of 2*pi/alat.
        """
        return self._kline

    @kline.setter
    def kline(self, kline):
        self._kline = numpy.array(kline)

        if len(self._kline.shape) != 1:
            raise ValueError(kline)

    @kline.deleter
    def kline(self):
        del self._kline

    @property
    def kweight(self):
        """A length-nkpoint ndarray of floats holding the weights of k-points,
    normalized to 1.
        """
        return self._kweight

    @kweight.setter
    def kweight(self, kweight):
        self._kweight = numpy.array(kweight, float)

        if len(self._kweight.shape) != 1:
            raise ValueError(kweight)

    @kweight.deleter
    def kweight(self):
        del self._kweight

    @property
    def nkpath(self):
        """An integer holding the number of vertices of the path in reciprocal
    space.
        """
        _list = []
        if hasattr(self, 'kpath'):
            _list.append(self._kpath.shape[0])
        if hasattr(self, 'kindex'):
            _list.append(self._kindex.shape[0])
        if hasattr(self, 'klabel'):
            _list.append(len(self._klabel))

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nkpath')

        return len(self._kpath)

    @property
    def kpath(self):
        """An nkpath by 3 ndarray of floats holding the coordinates of
    vertices of the path in reciprocal space, in units of kunit.
        """
        return self._kpath

    @kpath.setter
    def kpath(self, kpath):
        self._kpath = numpy.array(kpath, float)

        if len(self._kpath.shape) != 2 or self._kpath.shape[1] != 3:
            raise ValueError(kpath)

    @kpath.deleter
    def kpath(self):
        del self._kpath

    @property
    def kindex(self):
        """A length-nkpath ndarray of integers holding the indices of k-points
    along the path in reciprocal space corresponding to vertices.
        """
        return self._kindex

    @kindex.setter
    def kindex(self, kindex):
        self._kindex = numpy.array(kindex, int)

        if len(self._kindex.shape) != 1:
            raise ValueError(kindex)

    @kindex.deleter
    def kindex(self):
        del self._kindex

    @property
    def klabel(self):
        """A length-nkpath list of strings holding the labels of vertices of
    the path in reciprocal space.
        """
        return self._klabel

    @klabel.setter
    def klabel(self, klabel):
        self._klabel = klabel

        if type(self._klabel) is not list:
            raise ValueError(klabel)

    @klabel.deleter
    def klabel(self):
        del self._klabel

    def set_alat(self, alat):
        """Method for setting the new value of alat = float
    and recalculating properties of class Cell,
    kpoint, kline, and kpath.
        """
        if abs(alat - self.alat) > common.EPS12:
            if hasattr(self, 'kline'):
                self.kline *= alat / self.alat
            if self.kunit == 'cartesian':
                if hasattr(self, 'kpoint'):
                    self.kpoint *= alat / self.alat
                if hasattr(self, 'kpath'):
                    self.kpath *= alat / self.alat
        super().set_alat(alat)

    def set_kunit(self, kunit):
        """Method for setting the new value of kunit =
    'cartesian'|'crystal' and recalculating
    kpoint and kpath.
        """
        if kunit != self.kunit:
            if kunit == 'cartesian':
                if not hasattr(self, 'bvec'):
                    self.calc_bvec()
                mtrx = self.bvec
            else:
                if not hasattr(self, 'avec'):
                    self.calc_avec()
                mtrx = self.avec.transpose()
            if hasattr(self, 'kpoint'):
                self.kpoint = numpy.dot(self.kpoint, mtrx)
            if hasattr(self, 'kpath'):
                self.kpath = numpy.dot(self.kpath, mtrx)
            self.kunit = kunit

    def calc_kindex(self, mode):
        """Method for calculating numbers of k-points along
    different sections of the path in reciprocal space
    given the number of k-points along the first section,
    kindex[1] - kindex[0], where kindex[0] = 0.
    Other sections will have the same number (mode = 'number')
    or density (mode = 'density') of k-points.
        """
        nkpath = self.nkpath
        nkfirst = self.kindex[1]
        nksect = nkpath - 1
        kindex = numpy.empty((nkpath), int)

        if mode == 'number':
            for ii in range(2, nkpath):
                kindex[ii] = ii * nkfirst

        elif mode == 'density':
            oldkunit = self.kunit
            self.set_kunit('cartesian')
            klen = numpy.zeros(nksect, float)
            for ii in range(nksect):
                klen[ii] = numpy.linalg.norm(self.kpath[
                        ii + 1, :] - self.kpath[ii, :])
            self.set_kunit(oldkunit)
            knum = [nkfirst]
            for ii in range(1, nksect):
                knum.append(int(round(float(nkfirst) * klen[ii] / klen[0])))
            for ii in range(2, nkpath):
                kindex[ii] = sum(knum[:ii])

        else:
            raise ValueError(mode)

        self.kindex = kindex

    def check_kindex(self):
        """Method for checking numbers of k-points along
    different sections of the path in reciprocal space.
    Returned value = 'uninitialized'|'number'|'density',
    the meaning of 'number'|'density' is the same as in
    calc_kindex.
        """
        kindex = self.kindex
        nksect = self.nkpath - 1
        knum = []
        for ii in range(nksect):
            knum.append(kindex[ii + 1] - kindex[ii])

        if any(knum[ii] < 1 for ii in range(nksect)):
            mode = 'uninitialized'
        elif knum[1:] == knum[:-1]:
            mode = 'number'
        else:
            mode = 'density'

        return mode

    def calc_kpoint(self):
        """Method for calculating kpoint given kpath."""
        kpath = self.kpath
        kindex = self.kindex
        nksect = self.nkpath - 1
        nkpoint = kindex[nksect] + 1

        kpoint = numpy.empty((nkpoint, 3), float)
        for ii in range(nksect):
            kstart = kpath[ii]
            kend = kpath[ii + 1]
            nkstart = kindex[ii]
            nkend = kindex[ii + 1]
            kstep = (kend - kstart) / float(nkend - nkstart)
            kpoint[nkstart:nkend, :] = numpy.add(numpy.outer(
                    numpy.arange(nkend - nkstart), kstep), kstart)
        kpoint[nkpoint - 1, :] = kpath[nksect]

        self.kpoint = kpoint

    def calc_kline(self):
        """Method for calculating kline given kpoint."""
        oldkunit = self.kunit
        self.set_kunit('cartesian')
        kdiff = numpy.diff(self.kpoint, axis=0)
        self.set_kunit(oldkunit)
        knorm = numpy.linalg.norm(kdiff, axis=1)
        kcumsum = numpy.cumsum(knorm)
        self.kline = numpy.insert(kcumsum, 0, 0.0)

    def read(self, fileformat, filenames, lapwkunit=None):
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
    'wannier-in'     ['seedname.win']
    'wannier-out'    ['seedname.win', 'seedname_band.dat']
                  or ['seedname.win', 'seedname_band.datup',
                             'seedname_band.datdn']
    'vasp-out'       ['OUTCAR']
    'lapw-out'       ['case.output1']
                  or ['case.output1up', 'case.output1dn']

    lapwkunit = 'cartesian'|'crystal' used for 'lapw-kpt',
            WIEN2k requires k-points in crystal coordinates
            with respect to conventional reciprocal lattice
            vectors, see xcrysden/tests/supportInfo.kpath
            for details, use 'cartesian' for fcc or bcc and
            'crystal' for hcp
        """
        if lapwkunit is None:
            lapwkunit = 'cartesian'

        if fileformat == 'internal':
            self._read_file_internal(2, filenames)
        elif fileformat == 'pw-out':
            self._read_file_pw_out(2, filenames)
        elif fileformat == 'bands-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_bands_out(2, [filenames[1]])
        elif fileformat == 'matdyn-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_matdyn_out(2, [filenames[1]])
        elif fileformat == 'inteqp-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_inteqp_out(2, [filenames[1]], None)
        elif fileformat == 'sigma-out':
            self._read_file_pw_out(1, [filenames[0]])
            self._read_file_sigma_out(2, [filenames[1]], None)
        elif fileformat == 'wannier-in':
            self._read_file_wannier_in(2, filenames)
        elif fileformat == 'wannier-out':
            self._read_file_wannier_in(1, [filenames[0]])
            self._read_file_wannier_out(2, [filenames[1]])
        elif fileformat == 'vasp-out':
            self._read_file_vasp_out(2, filenames)
        elif fileformat == 'lapw-out':
            self._read_file_lapw_out(2, filenames, lapwkunit)
        else:
            raise ValueError(fileformat)

    def write(self, fileformat, filenames, lapwkunit=None):
        """Method for writing properties to file.

    fileformat       filenames
    ----------       ---------
    'internal'       ['prefix.brave']
    'pw-in'          ['prefix.in']
    'wannier-in'     ['seedname.win']
    'vasp-kpt'       ['KPOINTS']
    'lapw-kpt'       ['case.klist_band']

    lapwkunit = 'cartesian'|'crystal' used for 'lapw-kpt',
            WIEN2k requires k-points in crystal coordinates
            with respect to conventional reciprocal lattice
            vectors, see xcrysden/tests/supportInfo.kpath
            for details, use 'cartesian' for fcc or bcc and
            'crystal' for hcp
        """
        if lapwkunit is None:
            lapwkunit = 'cartesian'

        if fileformat == 'internal':
            self._write_file_internal(2, filenames)
        elif fileformat == 'pw-in':
            self._write_file_pw_in(2, filenames)
        elif fileformat == 'wannier-in':
            self._write_file_wannier_in(2, filenames)
        elif fileformat == 'vasp-kpt':
            self._write_file_vasp_kpt(2, filenames)
        elif fileformat == 'lapw-kpt':
            self._write_file_lapw_kpt(2, filenames, lapwkunit)
        else:
            raise ValueError(fileformat)

    def __init__(
            self, kunit=None, kpoint=None, kline=None, kweight=None,
            kpath=None, kindex=None, klabel=None, **kwargs):
        super().__init__(**kwargs)

        if kunit is not None:
            self.kunit = kunit
        if kpoint is not None:
            self.kpoint = kpoint
        if kline is not None:
            self.kline = kline
        if kweight is not None:
            self.kweight = kweight
        if kpath is not None:
            self.kpath = kpath
        if kindex is not None:
            self.kindex = kindex
        if klabel is not None:
            self.klabel = klabel

