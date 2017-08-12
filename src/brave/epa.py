"""This module defines class EPA."""

import math

import numpy

import brave.common as common
from brave.dos import DOS

class EPA(DOS):
    """Class for representing parameters of the
    EPA (electron-phonon-averaged) approximation.
    """

    @property
    def nepa(self):
        """Number of points for the grid of energy|mu|temp|invtau,
    int.
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
        """Number of energy windows, int."""
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
        """Number of phonon modes, int."""
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
        """Maximum number of energy bins in all energy windows, int."""
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
        """Values of electron energies, array of nepa floats,
    in units of eV.
        """
        return self._energy

    @energy.setter
    def energy(self, energy):
        self._energy = numpy.array(energy, float)

    @energy.deleter
    def energy(self):
        del self._energy

    @property
    def mu(self):
        """Values of chemical potentials of electrons, or Fermi
    levels, array of nepa floats, in units of eV.
        """
        return self._mu

    @mu.setter
    def mu(self, mu):
        self._mu = numpy.array(mu, float)

    @mu.deleter
    def mu(self):
        del self._mu

    @property
    def temp(self):
        """Values of temperatures, array of nepa floats,
    in units of K.
        """
        return self._temp

    @temp.setter
    def temp(self, temp):
        self._temp = numpy.array(temp, float)

    @temp.deleter
    def temp(self):
        del self._temp

    @property
    def ee(self):
        """Edges of the energy windows, array of nwin floats, eV."""
        return self._ee

    @ee.setter
    def ee(self, ee):
        self._ee = numpy.array(ee, float)

    @ee.deleter
    def ee(self):
        del self._ee

    @property
    def de(self):
        """Widths of the energy bins, array of nwin floats, eV."""
        return self._de

    @de.setter
    def de(self, de):
        self._de = numpy.array(de, float)

    @de.deleter
    def de(self):
        del self._de

    @property
    def ne(self):
        """Numbers of energy bins, array of nwin ints."""
        return self._ne

    @ne.setter
    def ne(self, ne):
        self._ne = numpy.array(ne, int)

    @ne.deleter
    def ne(self):
        del self._ne

    @property
    def wavg(self):
        """Averaged phonon frequencies, array of nmode floats,
    cm^-1.
        """
        return self._wavg

    @wavg.setter
    def wavg(self, wavg):
        self._wavg = numpy.array(wavg, float)

    @wavg.deleter
    def wavg(self):
        del self._wavg

    @property
    def gavg(self):
        """Averaged squared absolute electron-phonon coupling
    matrix elements, array of nmode by nemax by
    nemax by nwin floats, eV^2.
        """
        return self._gavg

    @gavg.setter
    def gavg(self, gavg):
        self._gavg = numpy.array(gavg, float)

    @gavg.deleter
    def gavg(self):
        del self._gavg

    @property
    def invtau(self):
        """Values of inverse electron relaxation times,
    array of nepa floats, in units of s^-1.
        """
        return self._invtau

    @invtau.setter
    def invtau(self, invtau):
        self._invtau = numpy.array(invtau, float)

    @invtau.deleter
    def invtau(self):
        del self._invtau

    def read(self, fileformat, filenames):
        """Method for reading properties from file.

    fileformat         filenames
    ----------         ---------
    'boltztrap-dos'    ['case.intrans', 'case.transdos']
    'matdyn-dos'       ['prefix.vdos']
    'epa-out'          ['epa.dat']

    Inherits fileformat and filenames from class DOS.
        """

        if fileformat.lower() == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames)
        elif fileformat.lower() == 'matdyn-dos':
            self._read_dos_matdyn_dos(filenames)
        elif fileformat.lower() == 'epa-out':
            self._read_epa_epa_out(filenames)
        else:
            super().read(fileformat, filenames)

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
                idx1, idx2, weight1, weight2 = common._int_pts(
                        dos[0], en[nn] + ww[ll])
                dosa = weight1 * dos[1, idx1] + weight2 * dos[1, idx2]
                idx1, idx2, weight1, weight2 = common._int_pts(
                        dos[0], en[nn] - ww[ll])
                dose = weight1 * dos[1, idx1] + weight2 * dos[1, idx2]

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

    def __init__(
            self, energy=None, mu=None, temp=None, ee=None, de=None, ne=None,
            wavg=None, gavg=None, invtau=None, **kwargs):
        super().__init__(**kwargs)

        if energy != None:
            self.energy = energy
        if mu != None:
            self.mu = mu
        if temp != None:
            self.temp = temp
        if ee != None:
            self.ee = ee
        if de != None:
            self.de = de
        if ne != None:
            self.ne = ne
        if wavg != None:
            self.wavg = wavg
        if gavg != None:
            self.gavg = gavg
        if invtau != None:
            self.invtau = invtau

    def _read_epa_epa_out(self, filenames):
        contents = common._read_file(filenames)

        nn = 0
        tt = contents[0][nn].split()
        nn += 1
        nwin = int(tt[0])
        nmode = int(tt[1])

        ee = numpy.empty(nwin, float)
        de = numpy.empty(nwin, float)
        ne = numpy.empty(nwin, int)
        for ii in range(nwin):
            tt = contents[0][nn].split()
            nn += 1
            ee[ii] = float(tt[0])
            de[ii] = float(tt[1])
            ne[ii] = int(tt[2])

        wavg = numpy.empty(nmode, float)
        tt = contents[0][nn].split()
        nn += 1
        for ii in range(nmode):
            wavg[ii] = float(tt[ii])

        nemax = numpy.amax(ne)
        gavg = numpy.zeros((nmode, nemax, nemax, nwin), float)
        for ii in range(nwin):
            for jj in range(ne[ii]):
                for kk in range(ne[ii]):
                    tt = contents[0][nn].split()
                    nn += 1
                    for ll in range(nmode):
                        gavg[ll, kk, jj, ii] = float(tt[ll + 3])

        self.ee, self.de, self.ne = ee, de, ne
        self.wavg, self.gavg = wavg, gavg

