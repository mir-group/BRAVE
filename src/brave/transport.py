"""This module defines class Transport."""

import numpy

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
    def mu(self, mu):
        self._mu = numpy.array(mu, float)

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
    def temp(self, temp):
        self._temp = numpy.array(temp, float)

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
    def numelec(self, numelec):
        self._numelec = numpy.array(numelec, float)

        if len(self._numelec.shape) != 2 or self._numelec.shape[
                0] != self.nmu or self._numelec.shape[1] != self.ntemp:
            raise ValueError(numelec)

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
    def convdos(self, convdos):
        self._convdos = numpy.array(convdos, float)

        if len(self._convdos.shape) != 2 or self._convdos.shape[
                0] != self.nmu or self._convdos.shape[1] != self.ntemp:
            raise ValueError(convdos)

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
    def seebeck(self, seebeck):
        self._seebeck = numpy.array(seebeck, float)

        if len(self._seebeck.shape) != 2 or self._seebeck.shape[
                0] != self.nmu or self._seebeck.shape[1] != self.ntemp:
            raise ValueError(seebeck)

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
    def sigma(self, sigma):
        self._sigma = numpy.array(sigma, float)

        if len(self._sigma.shape) != 2 or self._sigma.shape[
                0] != self.nmu or self._sigma.shape[1] != self.ntemp:
            raise ValueError(sigma)

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
    def hall(self, hall):
        self._hall = numpy.array(hall, float)

        if len(self._hall.shape) != 2 or self._hall.shape[
                0] != self.nmu or self._hall.shape[1] != self.ntemp:
            raise ValueError(hall)

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
    def kappael(self, kappael):
        self._kappael = numpy.array(kappael, float)

        if len(self._kappael.shape) != 2 or self._kappael.shape[
                0] != self.nmu or self._kappael.shape[1] != self.ntemp:
            raise ValueError(kappael)

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
    def specheat(self, specheat):
        self._specheat = numpy.array(specheat, float)

        if len(self._specheat.shape) != 2 or self._specheat.shape[
                0] != self.nmu or self._specheat.shape[1] != self.ntemp:
            raise ValueError(specheat)

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
    def magsus(self, magsus):
        self._magsus = numpy.array(magsus, float)

        if len(self._magsus.shape) != 2 or self._magsus.shape[
                0] != self.nmu or self._magsus.shape[1] != self.ntemp:
            raise ValueError(magsus)

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
    def kappalat(self, kappalat):
        self._kappalat = numpy.array(kappalat, float)

        if len(self._kappalat.shape) != 2 or self._kappalat.shape[
                0] != self.nmu or self._kappalat.shape[1] != self.ntemp:
            raise ValueError(kappalat)

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
    def kappa(self, kappa):
        self._kappa = numpy.array(kappa, float)

        if len(self._kappa.shape) != 2 or self._kappa.shape[
                0] != self.nmu or self._kappa.shape[1] != self.ntemp:
            raise ValueError(kappa)

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
    def L(self, L):
        self._L = numpy.array(L, float)

        if len(self._L.shape) != 2 or self._L.shape[
                0] != self.nmu or self._L.shape[1] != self.ntemp:
            raise ValueError(L)

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
    def PF(self, PF):
        self._PF = numpy.array(PF, float)

        if len(self._PF.shape) != 2 or self._PF.shape[
                0] != self.nmu or self._PF.shape[1] != self.ntemp:
            raise ValueError(PF)

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
    def ZT(self, ZT):
        self._ZT = numpy.array(ZT, float)

        if len(self._ZT.shape) != 2 or self._ZT.shape[
                0] != self.nmu or self._ZT.shape[1] != self.ntemp:
            raise ValueError(ZT)

    @ZT.deleter
    def ZT(self):
        del self._ZT

    def model_kappalat(self, kappalatvalue, tempvalue):
        """Method for modeling kappalat dependence on mu and temp.
    kappalat = kappalatvalue * tempvalue / temp if tempvalue > 0,
    kappalat = kappalatvalue if tempvalue <= 0, kappalatvalue
    in units of W/(m K), tempvalue in units of K.
        """
        kappalat = numpy.zeros((self.nmu, self.ntemp), float)

        for imu in range(self.nmu):
            for itemp in range(self.ntemp):
                kappalat[imu, itemp] = kappalatvalue
                if tempvalue > common.EPS12:
                    kappalat[imu, itemp] *= tempvalue / self.temp[itemp]

        self.kappalat = kappalat

    def calc_kappa(self):
        """Method for calculating kappa given kappael and kappalat."""

        self.kappa = numpy.add(self.kappael, self.kappalat)

    def calc_L(self):
        """Method for calculating L given sigma and kappael."""

        _sigma = numpy.maximum(self.sigma, common.EPS12)
        _temp = numpy.broadcast_to(self.temp, (self.nmu, self.ntemp))
        _sigma_temp = numpy.multiply(_sigma, _temp)
        self.L = numpy.divide(self.kappael, _sigma_temp)

    def calc_PF(self):
        """Method for calculating PF given sigma and seebeck."""

        _seebeck2 = numpy.square(self.seebeck)
        self.PF = numpy.multiply(self.sigma, _seebeck2)

    def calc_ZT(self):
        """Method for calculating ZT given PF and kappa."""

        _temp = numpy.broadcast_to(self.temp, (self.nmu, self.ntemp))
        _PF_temp = numpy.multiply(self.PF, _temp)
        self.ZT = numpy.divide(_PF_temp, self.kappa)

    def interpolate_binary(self, propname, paramtype, paramvalue):
        """Method for interpolating functional dependence of
    property propname on mu, temp, or numelec at a given
    value paramvalue = float of parameter paramtype.
    Returns array of nmu floats for paramtype = 'temp' or
    array of ntemp floats for paramtype = 'mu'|'numelec'.
    In case of paramtype = 'numelec' paramvalue can be
    array of ntemp floats.

    paramtype       paramvalue
    ---------       ----------
    'temp'          temperature in units of K
    'mu'            chemical potential in units of eV
    'numelec'       number of electrons in units of el/uc
                    (electrons per unit cell)
        """
        if propname.lower() == 'nmu' or propname.lower(
                ) == 'mu' or propname.lower() == 'ntemp' or propname.lower(
                ) == 'temp':
            raise ValueError(propname)

        propvalue = getattr(self, propname)

        if paramtype.lower() == 'temp':
            slice = numpy.zeros(self.nmu, float)
            itemp1, itemp2, weight1, weight2 = common._int_pts(
                    self.temp, paramvalue)
            for imu in range(self.nmu):
                slice[imu] = (weight1 * propvalue[
                        imu, itemp1] + weight2 * propvalue[imu, itemp2])
        elif paramtype.lower() == 'mu':
            slice = numpy.zeros(self.ntemp, float)
            imu1, imu2, weight1, weight2 = common._int_pts(
                    self.mu, paramvalue)
            for itemp in range(self.ntemp):
                slice[itemp] = (weight1 * propvalue[
                        imu1, itemp] + weight2 * propvalue[imu2, itemp])
        elif paramtype.lower() == 'numelec':
            slice = numpy.zeros(self.ntemp, float)
            numelec = self.numelec
            for itemp in range(self.ntemp):
                _paramvalue = paramvalue if type(
                        paramvalue) == float else paramvalue[itemp]
                slice[itemp] = numpy.interp(_paramvalue, numelec[
                        :, itemp], propvalue[:, itemp])
        else:
            raise ValueError(paramtype)

        return slice

    def interpolate_unary(self, propunary, argtype, argvalue):
        """Method for interpolating value of functional dependence
    of property propunary produced by method interpolate_binary
    at a given value argvalue = float of argument argtype.
    Returns float.

    argtype       argvalue
    -------       --------
    'temp'        temperature in units of K
    'mu'          chemical potential in units of eV
        """
        if argtype.lower() == 'temp':
            value = numpy.interp(argvalue, self.temp, propunary)
        elif argtype.lower() == 'mu':
            value = numpy.interp(argvalue, self.mu, propunary)
        else:
            raise ValueError(argtype)

        return value

    def renorm_numelec(self):
        """Method for renormalizing property numelec.
    Calculates numelec0, subtracts numelec0 from
    numelec, and returns numelec0. Here numelec0
    is interpolated from numelec, though it can
    be extracted from file case.transdos.
        """
        numelec = self.numelec
        numelec0 = numpy.interp(0.0, self.mu, numelec[:, 0])
        numelec -= numelec0
        self.numelec = numelec

        return numelec0

    def convert_numelec(self, inputvalue, inputunit, outputunit):
        """Method for converting number of electrons inputvalue =
    float from units of inputunit to units of outputunit.
    inputunit, outputunit = 'cm^-3'|'el/uc' (electrons
    per unit cell).
        """
        _unit_list = ['cm^-3', 'el/uc']
        if inputunit.lower() not in _unit_list:
            raise ValueError(inputunit)
        if outputunit.lower() not in _unit_list:
            raise ValueError(outputunit)
        outputvalue = inputvalue

        if outputunit.lower() != inputunit.lower():
            oldaunit = self.aunit
            self.set_aunit('angstrom')
            if outputunit.lower() == _unit_list[0]:
                outputvalue /= self.avol * self.alat ** 3 * 1.0e-24
            else:
                outputvalue *= self.avol * self.alat ** 3 * 1.0e-24
            self.set_aunit(oldaunit)

        return outputvalue

    def convert_argument(self, outputtype, inputvalues):
        """Method for converting between temp, mu, and numelec.

    outputtype       inputvalues
    ----------       -----------
    'temp'           [mu, numelec]
    'mu'             [temp, numelec]
    'numelec'        [temp, mu]

    temp in K, mu in eV, numelec in el/uc
    (electrons per unit cell)
        """
        if outputtype.lower() == 'temp':
            slice = self.interpolate_binary('numelec', 'mu', inputvalues[0])
            outputvalue = numpy.interp(inputvalues[1], slice, self.temp)
        elif outputtype.lower() == 'mu':
            slice = self.interpolate_binary('numelec', 'temp', inputvalues[0])
            outputvalue = numpy.interp(inputvalues[1], slice, self.mu)
        elif outputtype.lower() == 'numelec':
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
        """Method for optimizing mu within mu_interval at given
    temp to maximize property propname. Returns the optimal
    values of mu and numelec and the maximum value of
    propname. If fraction is present returns the intervals
    of mu and numelec for which propname is larger than the
    fraction of the maximum value of propname. mu_interval
    is a list of two floats in eV, temp is a float in K,
    fraction is a float between 0 and 1.
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
        if jmu == 0 or fraction == None:
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

        if fraction == None:
            return mu_opt, numelec_opt, propvalue_max
        else:
            return mu_opt, numelec_opt, propvalue_max, [mu_min, mu_max], [
                    numelec_min, numelec_max]

    def read(self, fileformat, filenames, tauvc=None, kappaelzeroj=False):
        """Method for reading properties from file.

    fileformat       filenames
    ----------       ---------
    'boltztrap-out'  ['case.intrans', 'case.trace']

    Inherits fileformat and filenames from class Cell.

    tauvc = electronic relaxation times for valence
    and conduction bands, list of two floats, in units
    of s, used for 'boltztrap-out'

    kappaelzeroj = True to compute kappael at zero electric
    current, kappael = kappa0 - sigma * seebeck^2 * temp,
    only for isotropic systems, used for 'boltztrap-out'
        """

        if fileformat.lower() == 'boltztrap-out':
            self._read_trn_boltztrap_out(filenames, tauvc, kappaelzeroj)
        else:
            super().read(fileformat, filenames)

    def __init__(
            self, mu=None, temp=None, numelec=None, convdos=None, seebeck=None,
            sigma=None, hall=None, kappael=None, specheat=None, magsus=None,
            kappalat=None, kappa=None, L=None, PF=None, ZT=None, **kwargs):
        super().__init__(**kwargs)

        if mu != None:
            self.mu = mu
        if temp != None:
            self.temp = temp
        if numelec != None:
            self.numelec = numelec
        if convdos != None:
            self.convdos = convdos
        if seebeck != None:
            self.seebeck = seebeck
        if sigma != None:
            self.sigma = sigma
        if hall != None:
            self.hall = hall
        if kappael != None:
            self.kappael = kappael
        if specheat != None:
            self.specheat = specheat
        if magsus != None:
            self.magsus = magsus
        if kappalat != None:
            self.kappalat = kappalat
        if kappa != None:
            self.kappa = kappa
        if L != None:
            self.L = L
        if PF != None:
            self.PF = PF
        if ZT != None:
            self.ZT = ZT

    def _read_trn_boltztrap_out(self, filenames, tauvc, kappaelzeroj):
        contents = common._read_file(filenames)

        efermi = float(contents[0][2].split()[0])

        nhead = 1
        ntot = len(contents[1]) - nhead
        for itemp in range(ntot):
            ftemp = float(contents[1][nhead + itemp].split()[1])
            if itemp > 0:
                if ftemp < ftemp_prev:
                    ntemp = itemp
                    break
            ftemp_prev = ftemp
        nmu = ntot // ntemp

        if hasattr(self, 'mu'):
            if nmu != self.nmu:
                raise ValueError('nmu mismatch')

        if hasattr(self, 'temp'):
            if ntemp != self.ntemp:
                raise ValueError('ntemp mismatch')

        mu = numpy.zeros(nmu, float)
        temp = numpy.zeros(ntemp, float)
        numelec = numpy.zeros((nmu, ntemp), float)
        convdos = numpy.zeros((nmu, ntemp), float)
        seebeck = numpy.zeros((nmu, ntemp), float)
        sigma = numpy.zeros((nmu, ntemp), float)
        hall = numpy.zeros((nmu, ntemp), float)
        kappael = numpy.zeros((nmu, ntemp), float)
        specheat = numpy.zeros((nmu, ntemp), float)
        magsus = numpy.zeros((nmu, ntemp), float)

        for imu in range(nmu):
            mu[imu] = (float(contents[1][nhead + imu * ntemp].split()[
                    0]) - efermi) * common.RYDBERG

        for itemp in range(ntemp):
            temp[itemp] = float(contents[1][nhead + itemp].split()[1])

        for imu in range(nmu):
            for itemp in range(ntemp):
                line = contents[1][nhead + imu * ntemp + itemp].split()
                numelec[imu, itemp] = -float(line[2])
                convdos[imu, itemp] = float(line[3]) / common.HARTREE
                seebeck[imu, itemp] = float(line[4])
                sigma[imu, itemp] = float(line[5])
                hall[imu, itemp] = float(line[6])
                kappael[imu, itemp] = float(line[7])
                specheat[imu, itemp] = float(line[8])
                magsus[imu, itemp] = float(line[9])

        if tauvc != None:
            for imu in range(nmu):
                if mu[imu] < 0.0:
                    tau = tauvc[0]
                else:
                    tau = tauvc[1]
                sigma[imu, :] *= tau
                kappael[imu, :] *= tau

        if kappaelzeroj:
            _seebeck2 = numpy.square(seebeck)
            _sigma_seebeck2 = numpy.multiply(sigma, _seebeck2)
            _temp = numpy.broadcast_to(temp, (nmu, ntemp))
            _sigma_seebeck2_temp = numpy.multiply(_sigma_seebeck2, _temp)
            kappael = numpy.subtract(kappael, _sigma_seebeck2_temp)

        self.mu = mu
        self.temp = temp
        self.numelec = numelec
        self.convdos = convdos
        self.seebeck = seebeck
        self.sigma = sigma
        self.hall = hall
        self.kappael = kappael
        self.specheat = specheat
        self.magsus = magsus

