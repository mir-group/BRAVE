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
    def ngrid(self):
        """An integer holding the number of energy grids."""
        _list = []
        if hasattr(self, 'edge'):
            _list.append(self._edge.shape[0])
        if hasattr(self, 'step'):
            _list.append(self._step.shape[0])
        if hasattr(self, 'nbin'):
            _list.append(self._nbin.shape[0])
        if hasattr(self, 'gavg'):
            _list.append(self._gavg.shape[0])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('ngrid')

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
    def nbinmax(self):
        """An integer holding the maximum number of energy bins in all energy
    grids.
        """
        _list = []
        if hasattr(self, 'gavg'):
            _list.append(self._gavg.shape[1])
            _list.append(self._gavg.shape[2])

        if _list[1:] == _list[:-1] and len(_list) > 0:
            return _list[0]
        else:
            raise AttributeError('nbinmax')

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
    def edge(self):
        """A length-ngrid ndarray of floats holding the grid edges of the
    energy grids, in units of eV.
        """
        return self._edge

    @edge.setter
    def edge(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('edge {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
            raise ValueError('edge {0!r}'.format(value))
        self._edge = value

    @edge.deleter
    def edge(self):
        del self._edge

    @property
    def step(self):
        """A length-ngrid ndarray of floats holding the grid steps of the
    energy grids, in units of eV.
        """
        return self._step

    @step.setter
    def step(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('step {0!r}'.format(value))
        if value.dtype != np.dtype('float') or len(value.shape) != 1:
            raise ValueError('step {0!r}'.format(value))
        self._step = value

    @step.deleter
    def step(self):
        del self._step

    @property
    def nbin(self):
        """A length-ngrid ndarray of integers holding the numbers of bins
    in the energy grids.
        """
        return self._nbin

    @nbin.setter
    def nbin(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('nbin {0!r}'.format(value))
        if value.dtype != np.dtype('int') or len(value.shape) != 1:
            raise ValueError('nbin {0!r}'.format(value))
        self._nbin = value

    @nbin.deleter
    def nbin(self):
        del self._nbin

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
        """A ngrid by nbinmax by nbinmax by nmode ndarray of floats holding
    the averaged squared absolute electron-phonon coupling matrix elements,
    in units of eV^2.
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

    Requires dunit and dos read with fileformat = 'boltztrap-dos' and edge,
    step, nbin, wavg and gavg read with fileformat = 'epa-out'.
        """
        if hasattr(self, 'invtau'):
            del self.invtau

        self.set_dunit(['uc', 'ev'])
        dos = self.dos

        nepa = self.nepa
        ngrid = self.ngrid
        nmode = self.nmode
        nbinmax = self.nbinmax
        en = self.energy
        mu = self.mu
        kt = self.temp * common.BOLTZMANN
        edge = self.edge
        step = self.step
        nbin = self.nbin
        # avoid the infrared divergence of the Bose-Einstein distribution
        ww = np.clip(self.wavg, 20.0, None) * 100 * common.PLANCK * common.LIGHT
        gavg = self.gavg

        gj = np.empty(nbinmax, float)
        gk = np.empty(2, float)
        invtau = np.empty(nepa, float)

        for nn in range(nepa):
            dummy = 0.0
            for ll in range(nmode):
                nw = 1 / (np.exp(ww[ll] / kt[nn]) - 1)
                fa = 1 / (np.exp((en[nn] + ww[ll] - mu[nn]) / kt[nn]) + 1)
                fe = 1 / (np.exp((en[nn] - ww[ll] - mu[nn]) / kt[nn]) + 1)
                dosa = np.interp(en[nn] + ww[ll], dos[0, :], dos[1, :])
                dose = np.interp(en[nn] - ww[ll], dos[0, :], dos[1, :])

                if (en[nn] < np.sum(edge) / ngrid):
                    ii = 0
                else:
                    ii = 1
                if (nbin[ii] == 1):
                    gk[:] = gavg[ii, 0, 0, ll]
                else:
                    xx = (en[nn] - edge[ii]) / step[ii]
                    xx = max(xx, common.EPS12)
                    xx = min(xx, nbin[ii] - common.EPS12)
                    jj = int(xx)
                    for kk in range(nbinmax):
                        gj[kk] = gavg[ii, jj, kk, ll]
                    for mm in range(2):
                        xx = (en[nn] + ww[ll] * (1 - 2 * mm) - edge[ii]) / step[
                            ii]
                        xx = max(xx, common.EPS12)
                        xx = min(xx, nbin[ii] - common.EPS12)
                        kk = int(xx)
                        gk[mm] = gj[kk]

                dummy += gk[0] * (nw + fa) * dosa + gk[1] * (nw + 1 - fe) * dose
            invtau[nn] = dummy

        self.invtau = invtau * math.pow(2.0 * math.pi, 2) / common.PLANCK

    def read(self, fileformat, filenames):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.

    fileformat         filenames
    ----------         ---------
    'epa-out'          ['case.intrans', 'prefix.epa.e']

    Inherits fileformat and filenames from class DOS.
        """

        if fileformat == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames)
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
            ngrid = int(tt[0])
            nmode = int(tt[1])

            edge = np.empty(ngrid, float)
            step = np.empty(ngrid, float)
            nbin = np.empty(ngrid, int)
            for ii in range(ngrid):
                tt = ff.readline().split()
                edge[ii] = float(tt[0])
                step[ii] = float(tt[1])
                nbin[ii] = int(tt[2])
            self.edge, self.step, self.nbin = (
                edge - efermi * common.RYDBERG, step, nbin)

            self.wavg = np.fromfile(
                ff, dtype = float, count = nmode, sep = ' ')

            nbinmax = np.amax(nbin)
            gavg = np.zeros((2, nbinmax, nbinmax, nmode), float)
            dummy = np.loadtxt(ff, dtype = float, usecols = (
                ii for ii in range(3, 3 + nmode)))
            gavg[0, :nbin[0], :nbin[0], :] = dummy[:nbin[0]**2, :].reshape(
                nbin[0], nbin[0], nmode)
            gavg[1, :nbin[1], :nbin[1], :] = dummy[nbin[0]**2:, :].reshape(
                nbin[1], nbin[1], nmode)
            self.gavg = gavg

    def __init__(
            self, energy = None, mu = None, temp = None, edge = None,
            step = None, nbin = None, wavg = None, gavg = None, invtau = None,
            **kwargs):
        super().__init__(**kwargs)

        if energy is not None:
            self.energy = energy
        if mu is not None:
            self.mu = mu
        if temp is not None:
            self.temp = temp
        if edge is not None:
            self.edge = edge
        if step is not None:
            self.step = step
        if nbin is not None:
            self.nbin = nbin
        if wavg is not None:
            self.wavg = wavg
        if gavg is not None:
            self.gavg = gavg
        if invtau is not None:
            self.invtau = invtau

