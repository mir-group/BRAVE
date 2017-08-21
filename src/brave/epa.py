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
    def energy(self, energy):
        self._energy = numpy.array(energy, float)

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
    def mu(self, mu):
        self._mu = numpy.array(mu, float)

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
    def temp(self, temp):
        self._temp = numpy.array(temp, float)

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
    def ee(self, ee):
        self._ee = numpy.array(ee, float)

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
    def de(self, de):
        self._de = numpy.array(de, float)

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
    def ne(self, ne):
        self._ne = numpy.array(ne, int)

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
    def wavg(self, wavg):
        self._wavg = numpy.array(wavg, float)

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
    def gavg(self, gavg):
        self._gavg = numpy.array(gavg, float)

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
    def invtau(self, invtau):
        self._invtau = numpy.array(invtau, float)

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
        if soc == None:
            soc = False

        if fileformat.lower() == 'boltztrap-dos':
            self._read_dos_boltztrap_dos(filenames, soc)
        elif fileformat.lower() == 'matdyn-dos':
            self._read_dos_matdyn_dos(filenames)
        elif fileformat.lower() == 'epa-out':
            self._read_epa_epa_out(filenames)
        else:
            super().read(fileformat, filenames)

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

