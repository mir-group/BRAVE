"""This module defines class EPA."""

import math

import numpy

import brave.common as common
from brave.dos import DOS

class EPA(DOS):
    """Class for representing the EPA electron relaxation time.

    Class EPA computes the electron relaxation time using the
    electron-phonon-averaged approximation.
    """

    @property
    def nepa(self):
        """An integer holding the number of grid points for energy, mu, temp
    and invtau.
        """
        _list = []
        if hasattr(self, 'energy'):
            _list.append(self._energy.shape[0])
        if hasattr(self, 'mu'):
            _list.append(self._mu.shape[0])
        if hasattr(self, 'temp'):
            _list.append(self._temp.shape[0])
        if hasattr(self, 'invtau'):
            _list.append(self._invtau.shape[0])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nepa')

    @property
    def nwin(self):
        """An integer holding the number of energy windows."""
        _list = []
        if hasattr(self, 'ee'):
            _list.append(self._ee.shape[0])
        if hasattr(self, 'de'):
            _list.append(self._de.shape[0])
        if hasattr(self, 'ne'):
            _list.append(self._ne.shape[0])
        if hasattr(self, 'gavg'):
            _list.append(self._gavg.shape[3])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nwin')

    @property
    def nmode(self):
        """An integer holding the number of phonon modes."""
        _list = []
        if hasattr(self, 'wavg'):
            _list.append(self._wavg.shape[0])
        if hasattr(self, 'gavg'):
            _list.append(self._gavg.shape[0])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nmode')

    @property
    def nemax(self):
        """An integer holding the maximum number of energy bins in all energy
    windows.
        """
        _list = []
        if hasattr(self, 'gavg'):
            _list.append(self._gavg.shape[1])
            _list.append(self._gavg.shape[2])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nemax')

    @property
    def energy(self):
        """A length-nepa ndarray of floats holding the electron energies, in
    units of eV.
        """
        return self._energy

    @energy.setter
    def energy(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('energy {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('energy {0!r}'.format(value))
        self._energy = value

    @energy.deleter
    def energy(self):
        del self._energy

    @property
    def mu(self):
        """A length-nepa ndarray of floats holding the chemical potentials of
    electrons or the Fermi levels, in units of eV.
        """
        return self._mu

    @mu.setter
    def mu(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('mu {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('mu {0!r}'.format(value))
        self._mu = value

    @mu.deleter
    def mu(self):
        del self._mu

    @property
    def temp(self):
        """A length-nepa ndarray of floats holding the temperatures, in units
    of K.
        """
        return self._temp

    @temp.setter
    def temp(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('temp {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('temp {0!r}'.format(value))
        self._temp = value

    @temp.deleter
    def temp(self):
        del self._temp

    @property
    def ee(self):
        """A length-nwin ndarray of floats holding the edges of the energy
    windows, in units of eV.
        """
        return self._ee

    @ee.setter
    def ee(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('ee {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('ee {0!r}'.format(value))
        self._ee = value

    @ee.deleter
    def ee(self):
        del self._ee

    @property
    def de(self):
        """A length-nwin ndarray of floats holding the widths of the energy
    bins, in units of eV.
        """
        return self._de

    @de.setter
    def de(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('de {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('de {0!r}'.format(value))
        self._de = value

    @de.deleter
    def de(self):
        del self._de

    @property
    def ne(self):
        """A length-nwin ndarray of integers holding the numbers of energy
    bins.
        """
        return self._ne

    @ne.setter
    def ne(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('ne {0!r}'.format(value))
        if value.dtype != numpy.dtype('int') or len(value.shape) != 1:
            raise ValueError('ne {0!r}'.format(value))
        self._ne = value

    @ne.deleter
    def ne(self):
        del self._ne

    @property
    def wavg(self):
        """A length-nmode ndarray of floats holding the averaged phonon
    frequencies, in units of cm^-1.
        """
        return self._wavg

    @wavg.setter
    def wavg(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('wavg {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('wavg {0!r}'.format(value))
        self._wavg = value

    @wavg.deleter
    def wavg(self):
        del self._wavg

    @property
    def gavg(self):
        """A nmode by nemax by nemax by nwin ndarray of floats holding the
    averaged squared absolute electron-phonon coupling matrix elements, in
    units of eV^2.
        """
        return self._gavg

    @gavg.setter
    def gavg(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('gavg {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 4:
            raise ValueError('gavg {0!r}'.format(value))
        self._gavg = value

    @gavg.deleter
    def gavg(self):
        del self._gavg

    @property
    def invtau(self):
        """A length-nepa ndarray of floats holding the inverse electron
    relaxation times, in units of s^-1.
        """
        return self._invtau

    @invtau.setter
    def invtau(self, value):
        if not isinstance(value, numpy.ndarray):
            raise TypeError('invtau {0!r}'.format(value))
        if value.dtype != numpy.dtype('float') or len(value.shape) != 1:
            raise ValueError('invtau {0!r}'.format(value))
        self._invtau = value

    @invtau.deleter
    def invtau(self):
        del self._invtau

    def calc_invtau(self):
        """Method for calculating invtau given energy, mu, temp.
    dunit, dos are read using fileformat = 'boltztrap-dos'.
    ee, de, ne, wavg, gavg are read using fileformat = 'epa-out'.
        """
        if hasattr(self, 'invtau'):
            del self.invtau

        self.set_dunit(['uc', 'ev'])
        dos = self.dos
        nspin = 2

        nepa = self.nepa
        nwin = self.nwin
        nmode = self.nmode
        nemax = self.nemax
        en = self.energy
        mu = self.mu
        kt = self.temp * common.BOLTZMANN
        ee = self.ee
        de = self.de
        ne = self.ne
        ww = self.wavg * 1.0e2 * common.PLANCK * common.LIGHT
        gavg = self.gavg

        gj = numpy.empty(nemax, float)
        gk = numpy.empty(2, float)
        invtau = numpy.empty(nepa, float)

        for nn in range(nepa):
            dummy = 0.0
            for ll in range(nmode):
                nw = 1.0 / (math.exp(ww[ll] / kt[nn]) - 1.0)
                fa = 1.0 / (math.exp((en[nn] + ww[ll] - mu[nn]) / kt[nn]
                        ) + 1.0)
                fe = 1.0 / (math.exp((en[nn] - ww[ll] - mu[nn]) / kt[nn]
                        ) + 1.0)
                dosa = numpy.interp(en[nn] + ww[ll], dos[0,:], dos[1,:])
                dose = numpy.interp(en[nn] - ww[ll], dos[0,:], dos[1,:])

                if (en[nn] < numpy.sum(ee) / nwin):
                    ii = 1
                else:
                    ii = 2
                if (ne[ii] == 1):
                    gk[:] = gavg(ll, 1, 1, ii)
                else:
                    xx = (en[nn] - ee[ii]) / de[ii]
                    xx = max(xx, 1.0e-12)
                    xx = min(xx, ne[ii] - 1.0e-12)
                    jj = int(xx)
                    for kk in range(nemax):
                        gj[kk] = gavg[ll, kk, jj, ii]
                    for mm in range(2):
                        xx = (en[nn] + ww[ll] * (3 - 2 * mm) - ee[ii]) / de[ii]
                        xx = max(xx, 1.0e-12)
                        xx = min(xx, ne[ii] - 1.0e-12)
                        kk = int(xx)
                        gk[mm] = gj[kk]

                dummy += gk[0] * (nw + fa) * dosa + gk[1] * (nw + 1.0 - fe
                        ) * dose
            invtau[nn] = dummy / nspin

        self.invtau = invtau * math.pow(2.0 * math.pi, 2) / common.PLANCK

    def read(self, fileformat, filenames, soc=None):
        """Method for reading properties from file.

    fileformat         filenames
    ----------         ---------
    'boltztrap-dos'    ['case.intrans', 'case.transdos']
    'matdyn-dos'       ['prefix.vdos']
    'epa-out'          ['epa.dat']

    Inherits fileformat and filenames from class DOS. Set soc
    to True if the calculation includes the spin-orbit coupling.
        """
        if soc is None:
            soc = False

        if fileformat == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames, soc)
        elif fileformat == 'matdyn-dos':
            self._read_dos_matdyn_dos(filenames)
        elif fileformat == 'epa-out':
            self._read_epa_epa_out(filenames)
        else:
            super().read(fileformat, filenames)

    def _read_epa_epa_out(self, filenames):

        with open(filenames[0]) as ff:
            tt = ff.readline().split()
            nwin = int(tt[0])
            nmode = int(tt[1])

            ee = numpy.empty(nwin, float)
            de = numpy.empty(nwin, float)
            ne = numpy.empty(nwin, int)
            for ii in range(nwin):
                tt = ff.readline().split()
                ee[ii] = float(tt[0])
                de[ii] = float(tt[1])
                ne[ii] = int(tt[2])

            wavg = numpy.empty(nmode, float)
            tt = ff.readline().split()
            for ii in range(nmode):
                wavg[ii] = float(tt[ii])

            nemax = numpy.amax(ne)
            gavg = numpy.zeros((nmode, nemax, nemax, nwin), float)
            for ii in range(nwin):
                for jj in range(ne[ii]):
                    for kk in range(ne[ii]):
                        tt = ff.readline().split()
                        for ll in range(nmode):
                            gavg[ll, kk, jj, ii] = float(tt[ll + 3])

        self.ee, self.de, self.ne, self.wavg, self.gavg = ee, de, ne, wavg, gavg

    def __init__(
            self, energy=None, mu=None, temp=None, ee=None, de=None, ne=None,
            wavg=None, gavg=None, invtau=None, **kwargs):
        super().__init__(**kwargs)

        if energy is not None:
            self.energy = energy
        if mu is not None:
            self.mu = mu
        if temp is not None:
            self.temp = temp
        if ee is not None:
            self.ee = ee
        if de is not None:
            self.de = de
        if ne is not None:
            self.ne = ne
        if wavg is not None:
            self.wavg = wavg
        if gavg is not None:
            self.gavg = gavg
        if invtau is not None:
            self.invtau = invtau

