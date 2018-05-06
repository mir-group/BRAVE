"""This module defines class Transport."""

import linecache
import numpy as np

import brave.common as common
from brave.cell import Cell

class Transport(Cell):
    """Class for representing the transport coefficients.

    Class Transport defines the electronic transport coefficients as a function
    of the chemical potential of electrons or the Fermi level and temperature.
    """

    @property
    def nmu(self):
        """An integer holding the number of grid points for mu.
        """
        return self._mu.shape[0]

    @property
    def mu(self):
        """A length-nmu ndarray of floats holding the chemical potentials of
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
    def ntemp(self):
        """An integer holding the number of grid points for temp."""
        return self._temp.shape[0]

    @property
    def temp(self):
        """A length-ntemp ndarray of floats holding the temperatures, in units
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
    def numelec(self):
        """A nmu by ntemp ndarray of floats holding the number of electrons
    relative to intrinsic material or the excess charge carrier concentration,
    in units of el/uc (electrons per unit cell). Positive for electron doping
    and negative for hole doping.
        """
        return self._numelec

    @numelec.setter
    def numelec(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('numelec {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('numelec {0!r}'.format(value))
        self._numelec = value

    @numelec.deleter
    def numelec(self):
        del self._numelec

    @property
    def convdos(self):
        """A nmu by ntemp ndarray of floats holding the convolution of DOS and
    [-df/de], in units of el/(uc eV) (electrons per unit cell per eV).
        """
        return self._convdos

    @convdos.setter
    def convdos(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('convdos {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('convdos {0!r}'.format(value))
        self._convdos = value

    @convdos.deleter
    def convdos(self):
        del self._convdos

    @property
    def seebeck(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the Seebeck tensor, in units of V/K.
        """
        return self._seebeck

    @seebeck.setter
    def seebeck(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('seebeck {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('seebeck {0!r}'.format(value))
        self._seebeck = value

    @seebeck.deleter
    def seebeck(self):
        del self._seebeck

    @property
    def sigma(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the electrical conductivity tensor, in units of 1/(ohm m).
        """
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('sigma {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('sigma {0!r}'.format(value))
        self._sigma = value

    @sigma.deleter
    def sigma(self):
        del self._sigma

    @property
    def hall(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the Hall tensor, in units of m^3/C.
        """
        return self._hall

    @hall.setter
    def hall(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('hall {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('hall {0!r}'.format(value))
        self._hall = value

    @hall.deleter
    def hall(self):
        del self._hall

    @property
    def kappael(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the electronic part of the thermal conductivity tensor, in units of
    W/(m K).
        """
        return self._kappael

    @kappael.setter
    def kappael(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('kappael {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('kappael {0!r}'.format(value))
        self._kappael = value

    @kappael.deleter
    def kappael(self):
        del self._kappael

    @property
    def specheat(self):
        """A nmu by ntemp ndarray of floats holding the electronic specific
    heat, in units of J/(mol K).
        """
        return self._specheat

    @specheat.setter
    def specheat(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('specheat {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('specheat {0!r}'.format(value))
        self._specheat = value

    @specheat.deleter
    def specheat(self):
        del self._specheat

    @property
    def magsus(self):
        """A nmu by ntemp ndarray of floats holding the Pauli magnetic
    susceptibility, in units of m^3/mol.
        """
        return self._magsus

    @magsus.setter
    def magsus(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('magsus {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('magsus {0!r}'.format(value))
        self._magsus = value

    @magsus.deleter
    def magsus(self):
        del self._magsus

    @property
    def kappalat(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the lattice part of the thermal conductivity tensor, in units of W/(m K).
        """
        return self._kappalat

    @kappalat.setter
    def kappalat(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('kappalat {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('kappalat {0!r}'.format(value))
        self._kappalat = value

    @kappalat.deleter
    def kappalat(self):
        del self._kappalat

    @property
    def kappa(self):
        """A nmu by ntemp ndarray of floats holding one third of the trace of
    the thermal conductivity tensor, in units of W/(m K).
        """
        return self._kappa

    @kappa.setter
    def kappa(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('kappa {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('kappa {0!r}'.format(value))
        self._kappa = value

    @kappa.deleter
    def kappa(self):
        del self._kappa

    @property
    def L(self):
        """A nmu by ntemp ndarray of floats holding the Lorenz number, in units
    of (W ohm)/K^2.
        """
        return self._L

    @L.setter
    def L(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('L {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('L {0!r}'.format(value))
        self._L = value

    @L.deleter
    def L(self):
        del self._L

    @property
    def PF(self):
        """A nmu by ntemp ndarray of floats holding the power factor, in units
    of W/(m K^2).
        """
        return self._PF

    @PF.setter
    def PF(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('PF {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('PF {0!r}'.format(value))
        self._PF = value

    @PF.deleter
    def PF(self):
        del self._PF

    @property
    def ZT(self):
        """A nmu by ntemp ndarray of floats holding the thermoelectric figure
    of merit ZT.
        """
        return self._ZT

    @ZT.setter
    def ZT(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError('ZT {0!r}'.format(value))
        if value.dtype != np.dtype('float') or value.shape != (
                self.nmu, self.ntemp):
            raise ValueError('ZT {0!r}'.format(value))
        self._ZT = value

    @ZT.deleter
    def ZT(self):
        del self._ZT

    def model_kappalat(self, kappalatvalue, tempvalue):
        """Models kappalat dependence on mu and temp.

    Args:
        kappalatvalue (float): Model parameter in units of W/(m K).
        tempvalue (float): Model parameter in units of K.

    if tempvalue >  0 then kappalat = kappalatvalue * tempvalue / temp,
    if tempvalue <= 0 then kappalat = kappalatvalue.
        """
        kappalat = np.zeros((self.nmu, self.ntemp), float)

        for imu in range(self.nmu):
            for itemp in range(self.ntemp):
                kappalat[imu, itemp] = kappalatvalue
                if tempvalue > common.EPS12:
                    kappalat[imu, itemp] *= tempvalue / self.temp[itemp]

        self.kappalat = kappalat

    def calc_kappa(self):
        """Calculates kappa from kappael and kappalat."""

        self.kappa = np.add(self.kappael, self.kappalat)

    def calc_L(self):
        """Calculates L from sigma and kappael."""

        _sigma = np.maximum(self.sigma, common.EPS12)
        _temp = np.broadcast_to(self.temp, (self.nmu, self.ntemp))
        _sigma_temp = np.multiply(_sigma, _temp)
        self.L = np.divide(self.kappael, _sigma_temp)

    def calc_PF(self):
        """Calculates PF from sigma and seebeck."""

        _seebeck2 = np.square(self.seebeck)
        self.PF = np.multiply(self.sigma, _seebeck2)

    def calc_ZT(self):
        """Calculates ZT from PF and kappa."""

        _temp = np.broadcast_to(self.temp, (self.nmu, self.ntemp))
        _PF_temp = np.multiply(self.PF, _temp)
        self.ZT = np.divide(_PF_temp, self.kappa)

    def interpolate_binary(self, propname, paramtype, paramvalue):
        """Interpolates property at some value of parameter.

    Args:
        propname (str): Name of property.
        paramtype (str): Name of parameter. Possible values are below.
        paramvalue (float or ndarray): Value of parameter. Possible values are
            below.

    Returns:
        propunary (ndarray): Interpolated property. Possible values are below.

    paramtype       paramvalue
    ---------       ----------
    'temp'          float, temperature in units of K
    'mu'            float, chemical potential in units of eV
    'numelec'       float or length-ntemp ndarray of floats, number of electrons
                        in units of el/uc (electrons per unit cell)

    paramtype       propunary
    ---------       ---------
    'temp'          length-nmu ndarray of floats
    'mu'            length-ntemp ndarray of floats
    'numelec'       length-ntemp ndarray of floats
        """
        if propname in ('nmu', 'mu', 'ntemp', 'temp'):
            raise ValueError(propname)

        propvalue = getattr(self, propname)

        if paramtype == 'temp':
            slice = np.zeros(self.nmu, float)
            itemp1, itemp2, weight1, weight2 = common._int_pts(
                    self.temp, paramvalue)
            for imu in range(self.nmu):
                slice[imu] = (weight1 * propvalue[
                        imu, itemp1] + weight2 * propvalue[imu, itemp2])
        elif paramtype == 'mu':
            slice = np.zeros(self.ntemp, float)
            imu1, imu2, weight1, weight2 = common._int_pts(
                    self.mu, paramvalue)
            for itemp in range(self.ntemp):
                slice[itemp] = (weight1 * propvalue[
                        imu1, itemp] + weight2 * propvalue[imu2, itemp])
        elif paramtype == 'numelec':
            slice = np.zeros(self.ntemp, float)
            numelec = self.numelec
            for itemp in range(self.ntemp):
                _paramvalue = paramvalue if isinstance(
                        paramvalue, float) else paramvalue[itemp]
                slice[itemp] = np.interp(_paramvalue, numelec[
                        :, itemp], propvalue[:, itemp])
        else:
            raise ValueError(paramtype)

        return slice

    def interpolate_unary(self, propunary, argtype, argvalue):
        """Interpolates the interpolated property at some value of argument.

    Args:
        propunary (ndarray): Interpolated property produced by method
            interpolate_binary.
        argtype (str): Name of argument. Possible values are below.
        argvalue (float): Value of argument. Possible values are below.

    Returns:
        propvalue (float): Interpolated property.

    argtype       argvalue
    -------       --------
    'temp'        temperature in units of K
    'mu'          chemical potential in units of eV
        """
        if argtype == 'temp':
            value = np.interp(argvalue, self.temp, propunary)
        elif argtype == 'mu':
            value = np.interp(argvalue, self.mu, propunary)
        else:
            raise ValueError(argtype)

        return value

    def renorm_numelec(self):
        """Renormalizes numelec by subtracting numelec0.

    Returns:
        numelec0 (float): Value of numelec at mu set to the middle of the band
            gap and temp = 0.

    Here numelec0 is interpolated from numelec, though it can be extracted from
    file case.transdos.
        """
        numelec = self.numelec
        numelec0 = np.interp(0.0, self.mu, numelec[:, 0])
        numelec -= numelec0
        self.numelec = numelec

        return numelec0

    def convert_numelec(self, inputvalue, inputunit, outputunit):
        """Converts the number of electrons between different units.

    Args:
        inputvalue (float): Number of electrons.
        inputunit (str): Units of inputvalue. Possible values are below.
        outputunit (str): Units of outputvalue. Possible values are below.

    Returns:
        outputvalue (float): Number of electrons.

    inputunit and outputunit are 'cm^-3' or 'el/uc' (electrons per unit cell).
        """
        _unit_list = ['cm^-3', 'el/uc']
        if inputunit not in _unit_list:
            raise ValueError(inputunit)
        if outputunit not in _unit_list:
            raise ValueError(outputunit)
        outputvalue = inputvalue

        if outputunit != inputunit:
            oldaunit = self.aunit
            self.set_aunit('angstrom')
            if outputunit == _unit_list[0]:
                outputvalue /= self.avol * self.alat ** 3 * 1.0e-24
            else:
                outputvalue *= self.avol * self.alat ** 3 * 1.0e-24
            self.set_aunit(oldaunit)

        return outputvalue

    def convert_argument(self, outputtype, inputvalues):
        """Converts between different arguments.

    Args:
        outputtype (str): Type of output parameter. Possible values are below.
        inputvalues (str): Values of input parameters. Possible values are
            below.

    Returns:
        outputvalue (float): Value of output parameter.

    outputtype       inputvalues
    ----------       -----------
    'temp'           [mu, numelec]
    'mu'             [temp, numelec]
    'numelec'        [temp, mu]

    temp is in units of K, mu is in units of eV, numelec is in units of el/uc
    (electrons per unit cell).
        """
        if outputtype == 'temp':
            slice = self.interpolate_binary('numelec', 'mu', inputvalues[0])
            outputvalue = np.interp(inputvalues[1], slice, self.temp)
        elif outputtype == 'mu':
            slice = self.interpolate_binary('numelec', 'temp', inputvalues[0])
            outputvalue = np.interp(inputvalues[1], slice, self.mu)
        elif outputtype == 'numelec':
            itemp1, itemp2, wtemp1, wtemp2 = common._int_pts(
                    self.temp, inputvalues[0])
            imu1, imu2, wmu1, wmu2 = common._int_pts(self.mu, inputvalues[1])
            outputvalue = wtemp1 * (wmu1 * self.numelec[
                    imu1, itemp1] + wmu2 * self.numelec[
                    imu2, itemp1]) + wtemp2 * (wmu1 * self.numelec[
                    imu1, itemp2] + wmu2 * self.numelec[imu2, itemp2])
        else:
            raise ValueError(outputtype)

        return outputvalue

    def optimize_mu(self, propname, mu_interval, temp, fraction=None):
        """Optimizes mu at a given temp to maximize some property.

    Args:
        propname (str): Name of property.
        mu_interval (list): Interval of mu values where to search for mu_opt.
        temp (float): Value of temperature.
        fraction (float): Fraction of propvalue_max for computing mu_range and
            numelec_range. Optional.

    Returns:
        mu_opt (float): Optimal value of mu in units of eV.
        numelec_opt (float): Optimal value of numelec.
        propvalue_max (float): Maximum value of propname.
        mu_range (list): Range of mu values for which propname is larger than
            propvalue_max times fraction in units of eV. Optional.
        numelec_range (list): Range of numelec values for which propname is
            larger than propvalue_max times fraction. Optional.

    mu_interval, mu_opt and mu_range are in units of eV, temp is in units of K,
    numelec_opt and numelec_range are in units of el/uc (electrons per unit
    cell).
        """
        nmu = self.nmu
        mu = self.mu
        numelec = self.interpolate_binary('numelec', 'temp', temp)
        propvalue = self.interpolate_binary(propname, 'temp', temp)

        jmu = 0
        propvalue_max = 0.0
        for imu in range(1, nmu - 1):
            if mu[imu] < mu_interval[0] or mu[imu] > mu_interval[1]:
                continue
            if propvalue[imu] >= propvalue[imu - 1] and propvalue[imu
                    ] >= propvalue[imu + 1] and propvalue[imu] > propvalue_max:
                propvalue_max = propvalue[imu]
                jmu = imu

        # refine propvalue_max by parabolic fit of propvalue[imu]
        #     around imu=jmu
        # assumes mu[imu] is equally spaced
        # yy = aa xx^2 + bb xx + cc
        # xx = -1, 0, 1 for mu[jmu - 1], mu[jmu], mu[jmu + 1]
        # skip interpolation if we are not near the maximum (if |x| > 1),
        #     for example propvalue_max is found at the edge of mu_interval
        if jmu == 0:
            mu_opt = 0.0
            numelec_opt = 0.0
            propvalue_max = 0.0
        else:
            cc = propvalue[jmu]
            bb = 0.5 * (propvalue[jmu + 1] - propvalue[jmu - 1])
            aa = 0.5 * (propvalue[jmu + 1] + propvalue[jmu - 1]) - cc
            xx = -0.5 * bb / aa
            if abs(xx) > 1.0:
                xx = 0.0
            mu_opt = mu[jmu] + (mu[jmu + 1] - mu[jmu]) * xx
            numelec_opt = numelec[jmu] + (numelec[jmu + 1] - numelec[jmu]) * xx
            propvalue_max = aa * xx ** 2 + bb * xx + cc

        # find intervals by linear fit of propvalue[imu] around imu=kmu
        if jmu == 0 or fraction is None:
            mu_min = mu_opt
            mu_max = mu_opt
            numelec_min = numelec_opt
            numelec_max = numelec_opt
        else:
            for imu in range(jmu, -1, -1):
                if propvalue[imu] < fraction * propvalue_max:
                    kmu = imu
                    break
            xx = (fraction * propvalue_max - propvalue[kmu]) / (
                    propvalue[kmu + 1] - propvalue[kmu])
            mu_min = mu[kmu] + (mu[kmu + 1] - mu[kmu]) * xx
            numelec_min = numelec[kmu] + (
                    numelec[kmu + 1] - numelec[kmu]) * xx

            for imu in range(jmu, nmu, 1):
                if propvalue[imu] < fraction * propvalue_max:
                    kmu = imu
                    break
            xx = (fraction * propvalue_max - propvalue[kmu - 1]) / (
                    propvalue[kmu] - propvalue[kmu - 1])
            mu_max = mu[kmu - 1] + (mu[kmu] - mu[kmu - 1]) * xx
            numelec_max = numelec[kmu - 1] + (
                    numelec[kmu] - numelec[kmu - 1]) * xx

        if fraction is None:
            return mu_opt, numelec_opt, propvalue_max
        else:
            return mu_opt, numelec_opt, propvalue_max, [mu_min, mu_max], [
                    numelec_min, numelec_max]

    def read(self, fileformat, filenames, tauvc=None, kappaelzeroj=None):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.
        tauvc (list): Constant electronic relaxation times for valence and
            conduction bands in units of s. Optional.
        kappaelzeroj (bool): Set to True to compute kappael at zero electric
            current, kappael = kappa0 - sigma * seebeck^2 * temp, only valid
            for isotropic systems. Optional.

    fileformat       filenames
    ----------       ---------
    'boltztrap-out'  ['case.intrans', 'case.trace']

    Inherits fileformat and filenames from class Cell.
        """
        if kappaelzeroj is None:
            kappaelzeroj = False

        if fileformat == 'boltztrap-out':
            self._read_trn_boltztrap_out(filenames, tauvc, kappaelzeroj)
        else:
            super().read(fileformat, filenames)

    def _read_trn_boltztrap_out(self, filenames, tauvc, kappaelzeroj):

        efermi = float(linecache.getline(filenames[0], 3).split()[0])
        trace = np.loadtxt(filenames[
            1], dtype = float, skiprows = 1, unpack = True)

        ntemp = np.where(trace[1, :] == trace[1, 0])[0][1]
        nmu = trace.shape[1] // ntemp

        self.mu = (trace[0, ::ntemp] - efermi) * common.RYDBERG
        self.temp = trace[1, 0:ntemp]
        view = trace.reshape(trace.shape[0], nmu, ntemp)
        self.numelec = -view[2, :, :]
        self.convdos = view[3, :, :] / common.HARTREE
        self.seebeck = view[4, :, :]
        self.sigma = view[5, :, :]
        self.hall = view[6, :, :]
        self.kappael = view[7, :, :]
        self.specheat = view[8, :, :]
        self.magsus = view[9, :, :]

        if tauvc is not None:
            tau = np.broadcast_to(np.where(self.mu < 0.0, tauvc[0], tauvc[
                1]).reshape(nmu, 1), (nmu, ntemp))
            self.sigma = np.multiply(self.sigma, tau)
            self.kappael = np.multiply(self.kappael, tau)

        if kappaelzeroj:
            self.kappael = np.subtract(self.kappael, np.multiply(
                np.multiply(self.sigma, np.square(
                self.seebeck)), np.broadcast_to(self.temp, (nmu, ntemp))))

    def __init__(
            self, mu=None, temp=None, numelec=None, convdos=None, seebeck=None,
            sigma=None, hall=None, kappael=None, specheat=None, magsus=None,
            kappalat=None, kappa=None, L=None, PF=None, ZT=None, **kwargs):
        super().__init__(**kwargs)

        if mu is not None:
            self.mu = mu
        if temp is not None:
            self.temp = temp
        if numelec is not None:
            self.numelec = numelec
        if convdos is not None:
            self.convdos = convdos
        if seebeck is not None:
            self.seebeck = seebeck
        if sigma is not None:
            self.sigma = sigma
        if hall is not None:
            self.hall = hall
        if kappael is not None:
            self.kappael = kappael
        if specheat is not None:
            self.specheat = specheat
        if magsus is not None:
            self.magsus = magsus
        if kappalat is not None:
            self.kappalat = kappalat
        if kappa is not None:
            self.kappa = kappa
        if L is not None:
            self.L = L
        if PF is not None:
            self.PF = PF
        if ZT is not None:
            self.ZT = ZT

