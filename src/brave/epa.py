"""This module defines class EPA."""

import math
import linecache

import numpy as np

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
            _list.append(self._gavg.shape[0])

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
            _list.append(self._gavg.shape[3])

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
        if not isinstance(value, np.ndarray):
            raise TypeError('energy {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('mu {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('temp {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('ee {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('de {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('ne {0!r}'.format(value))
        if value.dtype != np.dtype('int') or len(value.shape) != 1:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('wavg {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
            raise ValueError('wavg {0!r}'.format(value))
        self._wavg = value

    @wavg.deleter
    def wavg(self):
        del self._wavg

    @property
    def gavg(self):
        """A nwin by nemax by nemax by nmode ndarray of floats holding the
    averaged squared absolute electron-phonon coupling matrix elements, in
    units of eV^2.
        """
        return self._gavg

    @gavg.setter
    def gavg(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('gavg {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 4:
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
        if not isinstance(value, np.ndarray):
            raise TypeError('invtau {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
            raise ValueError('invtau {0!r}'.format(value))
        self._invtau = value

    @invtau.deleter
    def invtau(self):
        del self._invtau

    def calc_invtau(self):
        """Calculates invtau for energy, mu and temp.

    Requires dunit and dos read with fileformat = 'boltztrap-dos' and ee, de,
    ne, wavg and gavg read with fileformat = 'epa-out'.
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
        ww = self.wavg * 100 * common.PLANCK * common.LIGHT
        gavg = self.gavg

        gj = np.empty(nemax, float)
        gk = np.empty(2, float)
        invtau = np.empty(nepa, float)

        for nn in range(nepa):
            dummy = 0.0
            for ll in range(nmode):
                nw = 1 / (np.exp(ww[ll] / kt[nn]) - 1)
                fa = 1 / (np.exp((en[nn] + ww[ll] - mu[nn]) / kt[nn]) + 1)
                fe = 1 / (np.exp((en[nn] - ww[ll] - mu[nn]) / kt[nn]) + 1)
                dosa = np.interp(en[nn] + ww[ll], dos[0,:], dos[1,:])
                dose = np.interp(en[nn] - ww[ll], dos[0,:], dos[1,:])

                if (en[nn] < np.sum(ee) / nwin):
                    ii = 0
                else:
                    ii = 1
                if (ne[ii] == 1):
                    gk[:] = gavg[ii, 0, 0, ll]
                else:
                    xx = (en[nn] - ee[ii]) / de[ii]
                    xx = max(xx, common.EPS12)
                    xx = min(xx, ne[ii] - common.EPS12)
                    jj = int(xx)
                    for kk in range(nemax):
                        gj[kk] = gavg[ii, jj, kk, ll]
                    for mm in range(2):
                        xx = (en[nn] + ww[ll] * (1 - 2 * mm) - ee[ii]) / de[ii]
                        xx = max(xx, common.EPS12)
                        xx = min(xx, ne[ii] - common.EPS12)
                        kk = int(xx)
                        gk[mm] = gj[kk]

                dummy += gk[0] * (nw + fa) * dosa + gk[1] * (nw + 1 - fe) * dose
            invtau[nn] = dummy / nspin

        self.invtau = invtau * math.pow(2.0 * math.pi, 2) / common.PLANCK

    def read(self, fileformat, filenames, soc=None):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.
        soc (bool): Set to True if the calculation includes the spin-orbit
            coupling.

    fileformat         filenames
    ----------         ---------
    'epa-out'          ['case.intrans', 'epa.dat']

    Inherits fileformat and filenames from class DOS.
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

        efermi = float(linecache.getline(filenames[0], 3).split()[0])

        with open(filenames[1], 'rb') as ff:
            tt = ff.readline().split()
            nwin = int(tt[0])
            nmode = int(tt[1])

            ee = np.empty(nwin, float)
            de = np.empty(nwin, float)
            ne = np.empty(nwin, int)
            for ii in range(nwin):
                tt = ff.readline().split()
                ee[ii] = float(tt[0])
                de[ii] = float(tt[1])
                ne[ii] = int(tt[2])
            self.ee, self.de, self.ne = ee - efermi * common.RYDBERG, de, ne

            self.wavg = np.fromfile(
                ff, dtype = float, count = nmode, sep = ' ')

            nemax = np.amax(ne)
            gavg = np.zeros((2, nemax, nemax, nmode), float)
            dummy = np.loadtxt(ff, dtype = float, usecols = (
                ii for ii in range(3, 3 + nmode)))
            gavg[0, :ne[0], :ne[0], :] = dummy[:ne[0]**2, :].reshape(
                ne[0], ne[0], nmode)
            gavg[1, :ne[1], :ne[1], :] = dummy[ne[0]**2:, :].reshape(
                ne[1], ne[1], nmode)
            self.gavg = gavg

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

