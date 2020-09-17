"""This module defines class File."""

import sys
import math
import fractions
import io
import os

import numpy as np

import brave.common as common

class File(object):
    """Class providing common I/O methods.

    Class File contains common input and output methods for classes Cell,
    Kpoint, and Energy.

    Args:
        level (int): Indicates whether to read/write the attributes of each
            class.

    level > 0: read/write the attributes of class Cell.
    level > 1: read/write the attributes of class Kpoint.
    level > 2: read/write the attributes of class Energy.
    """

    def _read_file_internal(self, level, filenames):

        with open(filenames[0], 'rb') as ff:
            for line in ff:
                if line.strip().startswith(b'#'):
                    continue

                if level > 0:
                    if b'prefix' in line:
                        self.prefix = line[line.find(b'=') + 1:].strip().decode(
                            )
                    elif b'aunit' in line:
                        self.aunit = line[line.find(b'=') + 1:].strip().decode()
                    elif b'alat' in line:
                        self.alat = float(line[line.find(b'=') + 1:].strip())
                    elif b'avec' in line:
                        self.avec = np.genfromtxt(
                            ff, dtype = float, max_rows = 3)
                    elif b'bvec' in line:
                        self.bvec = np.genfromtxt(
                            ff, dtype = float, max_rows = 3)
                    elif b'avol' in line:
                        self.avol = float(line[line.find(b'=') + 1:].strip())
                    elif b'bvol' in line:
                        self.bvol = float(line[line.find(b'=') + 1:].strip())
                    elif b'natom' in line:
                        self.natom = int(line[line.find(b'=') + 1:].strip())
                    elif b'nelec' in line:
                        self.nelec = float(line[line.find(b'=') + 1:].strip())
                    elif b'sym' in line:
                        nsym = int(line.split()[1])
                        self.rot = np.genfromtxt(
                            ff, dtype = int, max_rows = nsym).reshape(
                            nsym, 3, 3)

                if level > 1:
                    if b'kunit' in line:
                        self.kunit = line[line.find(b'=') + 1:].strip().decode()
                    elif b'kpoint' in line:
                        nkpoint = int(line.split()[1])
                        self.kpoint = np.genfromtxt(
                            ff, dtype = float, max_rows = nkpoint)
                    elif b'kline' in line:
                        nkpoint = int(line.split()[1])
                        self.kline = np.fromfile(
                            ff, dtype = float, count = nkpoint, sep = ' ')
                    elif b'kweight' in line:
                        nkpoint = int(line.split()[1])
                        self.kweight = np.fromfile(
                            ff, dtype = float, count = nkpoint, sep = ' ')
                    elif b'kpath' in line:
                        nkpath = int(line.split()[1])
                        self.kpath = np.genfromtxt(
                            ff, dtype = float, max_rows = nkpath)
                    elif b'kindex' in line:
                        nkpath = int(line.split()[1])
                        self.kindex = np.fromfile(
                            ff, dtype = int, count = nkpath, sep = ' ')
                    elif b'klabel' in line:
                        nkpath = int(line.split()[1])
                        self.klabel = list(np.genfromtxt(
                            ff, dtype = str, max_rows = nkpath))

                if level > 2:
                    if b'eunit' in line:
                        self.eunit = line[line.find(b'=') + 1:].strip().decode()
                    elif b'energy' in line:
                        words = line.split()
                        nkpoint = int(words[1])
                        nband = int(words[2])
                        nspin = int(words[3])
                        self.energy = np.fromfile(
                            ff, dtype = float, count = (
                            nkpoint * nband * nspin), sep = ' ').reshape(
                            nkpoint, nband, nspin)
                    elif b'efermi' in line:
                        self.efermi = float(line[line.find(b'=') + 1:].strip())
                    elif b'vref' in line:
                        self.vref = float(line[line.find(b'=') + 1:].strip())

    def _read_file_pw_out(self, level, filenames):
        if level > 0:
            buf = io.BytesIO()
        if level > 1:
            nspin = 1

        with open(filenames[0], 'rb') as ff:
            ii = 0
            for line in ff:
                if level > 0:
                    if b'Writing output data file' in line:
                        ss = line.split()[4].decode()
                        self.prefix = ss[:ss.find('.')]
                    elif b'lattice parameter' in line:
                        self.aunit = 'bohr'
                        self.alat = float(line.split()[4])
                    elif b'crystal axes' in line:
                        self.avec = np.genfromtxt(
                            ff, dtype = float, usecols = (
                            3, 4, 5), max_rows = 3)
                    elif b'reciprocal axes' in line:
                        self.bvec = np.genfromtxt(
                            ff, dtype = float, usecols = (
                            3, 4, 5), max_rows = 3)
                    elif b'unit-cell volume' in line:
                        self.avol = float(line.split()[3]) / self.alat ** 3
                    elif b'number of atoms/cell' in line:
                        self.natom = int(line.split()[4])
                    elif b'number of electrons' in line:
                        self.nelec = float(line.split()[4])
                    elif b'Sym. Ops.' in line or b'Sym.Ops.' in line:
                        nsym = int(line.split()[0])
                    elif b'No symmetry found' in line:
                        nsym = 1
                    elif b'cryst.   s' in line:
                        for kk in range(3):
                            buf.write(line[19:53])
                            line = ff.readline()

                if level > 1:
                    if b'magnetic structure' in line:
                        nspin = 2
                    elif b'number of k points=' in line:
                        self.kunit = 'crystal'
                        nkpoint = int(line.split()[4])
                    elif b' cryst. coord.' in line:
                        kk = np.genfromtxt(ff, dtype = float, delimiter = (
                            20, 12, 12, 12, 7, 12), usecols = (
                            1, 2, 3, 5), max_rows = nkpoint)
                        self.kpoint, self.kweight = kk[:, :-1], kk[:, 3] / 2

                if level > 2:
                    if b'number of Kohn-Sham states' in line:
                        self.eunit = 'ev'
                        nband = int(line.split()[4])
                    elif b'FFT' in line:
                        energy = np.empty((nkpoint, nband, nspin), float)
                    elif b'bands (ev)' in line or b'band energies (ev)' in line:
                        energy[ii % nkpoint, :, ii // nkpoint] = np.genfromtxt(
                            ff, dtype = float, max_rows = nband // 8 + 1,
                            skip_header = 1, delimiter = (
                            11, 9, 9, 9, 9, 9, 9, 9)).flatten()[:nband]
                        ii += 1
                    elif b'the Fermi energy is' in line:
                        self.efermi = float(line.split()[4])
                    elif b'highest occupied, lowest unoccupied level' in line:
                        words = line.split()
                        self.efermi = (float(words[6]) + float(words[7])) / 2.0

            if level > 0:
                buf.seek(0)
                self.rot = np.loadtxt(buf, dtype = int).reshape(nsym, 3, 3)
                buf.close()
            if level > 2:
                self.energy = energy

    def _read_file_wannier_in(self, level, filenames):
        if level > 1:
            nkpersect = 100

        with open(filenames[0], 'rb') as ff:
            for line in ff:
                if level > 0:
                    if b'Begin Unit_Cell_Cart' in line:
                        self.aunit = ff.readline().strip().lower().decode()
                        self.avec = np.genfromtxt(
                            ff, dtype = float, max_rows = 3)
                        self.alat = 1.0

                if level > 1:
                    if b'bands_num_points' in line:
                        nkpersect = int(line.split()[2])
                    elif b'kpoint_path' in line:
                        buf = io.BytesIO()
                        for line in ff:
                            if b'kpoint_path' in line:
                                break
                            buf.write(line)
                        buf.seek(0)
                        crd = np.loadtxt(buf, dtype = float, usecols = (
                            0, 1, 2, 4, 5, 6))
                        buf.seek(0)
                        lbl = np.loadtxt(buf, dtype = 'S20', usecols = (3, 7))
                        buf.close()

                        if lbl[1:, 0] != lbl[:-1, 1] or crd[1:, 0:3] != crd[
                                :-1, 3:6]:
                            raise ValueError(kpath)

                        self.kpath = np.concatenate((crd[:, 0:3], crd[
                            -1:, 3:6]))
                        self.klabel = np.concatenate((lbl[:, 0], lbl[
                            -1:, 1])).tolist()
                        self.kindex = [nkpersect] + [0 for ii in range(
                            crd.shape[0])]
                        self.kunit = 'crystal'

    def _read_file_vasp_out(self, level, filenames):
        if level > 2:
            nspin, kk, ll = 1, 0, 0

        with open(filenames[0], 'rb') as ff:
            for line in ff:
                if level > 0:
                    if b'SYSTEM =' in line:
                        self.prefix = line[line.find(b'=') + 1:].strip().decode(
                            )
                    elif b'ALAT' in line:
                        self.alat = float(line[line.find(b'=') + 1:].strip())
                        self.aunit = 'angstrom'
                    elif b'Lattice vectors:' in line:
                        line = ff.readline()
                        buf = io.BytesIO()
                        for jj in range(3):
                            line = ff.readline().replace(b',', b'')
                            buf.write(line[line.find(b'(') + 1:line.find(b')')])
                        buf.seek(0)
                        self.avec = np.loadtxt(buf, dtype = float).reshape(
                                3, 3) / self.alat
                        buf.close()
                    elif b'NELECT' in line:
                        self.nelec = float(line.split()[2])

                if level > 1:
                    if b'Dimension of arrays:' in line:
                        line = ff.readline()
                        words = line.split()
                        nkpoint = int(words[words.index(b'NKPTS') + 2])
                        nband = int(words[words.index(b'NBANDS=') + 1])
                    elif b'k-points in reciprocal lattice and weights' in line:
                        dummy = np.genfromtxt(
                            ff, dtype = float, max_rows = nkpoint)
                        self.kpoint, self.kweight = dummy[:, :-1], dummy[:, 3]
                        self.kunit = 'crystal'

                if level > 2:
                    if b'ISPIN' in line:
                        nspin = int(line.split()[2])
                    elif b'E-fermi :' in line:
                        self.efermi = float(line.split()[2])
                    elif b'band No.  band energies     occupation' in line:
                        if kk == 0 and ll == 0:
                            energy = np.empty((nkpoint, nband, nspin), float)
                        energy[kk, :, ll] = np.genfromtxt(
                            ff, dtype = float, usecols = 1, max_rows = nband)
                        kk += 1
                        if kk == nkpoint:
                            kk = 0
                            ll += 1
                        if ll == nspin:
                            self.energy, self.eunit = energy, 'ev'

    def _calc_lapw_avec(self, lattype, apar, bpar, cpar, alpha, beta, gamma):

        if lattype[0:1] == b'P':
            phi = math.acos((math.cos(gamma) - math.cos(alpha) * math.cos(
                    beta)) / (math.sin(alpha) * math.sin(beta)))
            avec = np.array([[math.sin(beta) * math.sin(phi) * apar, math.sin(
                beta) * math.cos(phi) * apar, math.cos(beta) * apar], [
                0, math.sin(alpha) * bpar, math.cos(alpha) * bpar], [
                0, 0, cpar]], float)

        elif lattype[0:1] == b'H':
            avec = np.array([[0.5 * math.sqrt(3) * apar, -0.5 * apar, 0], [
                0, apar, 0], [0, 0, cpar]], float)

        elif lattype[0:1] in [b'T', b'R']:
            avec = np.array([[apar / (2 * math.sqrt(3)), -apar / 2, cpar / 3], [
                apar / (2 * math.sqrt(3)), apar / 2, cpar / 3], [
                -apar / math.sqrt(3), 0, cpar / 3]], float)

        elif lattype[0:1] == b'F':
            avec = np.array([[0, bpar / 2, cpar / 2], [apar / 2, 0, cpar / 2], [
                apar / 2, bpar / 2, 0]], float)

        elif lattype[0:1] == b'B':
            avec = np.array([[-apar / 2, bpar / 2, cpar / 2], [
                apar / 2, -bpar / 2, cpar / 2], [
                apar / 2, bpar / 2, -cpar / 2]], float)

        elif lattype[0:3] == b'CXY' and abs(
                gamma - math.pi / 2) < common.EPS12 and (abs(
                alpha - math.pi / 2) < common.EPS12 or abs(
                beta - math.pi / 2) < common.EPS12):
            avec = np.array([[math.sin(beta) * apar / 2, -math.sin(
                alpha) * bpar / 2, math.cos(beta) * apar / 2 - math.cos(
                alpha) * bpar / 2], [math.sin(beta) * apar / 2, math.sin(
                alpha) * bpar / 2, math.cos(beta) * apar / 2 + math.cos(
                alpha) * bpar / 2], [0, 0, cpar]], float)

        elif lattype[0:3] == b'CYZ' and abs(
                alpha - math.pi / 2) < common.EPS12 and (abs(
                beta - math.pi / 2) < common.EPS12 or abs(
                gamma - math.pi / 2) < common.EPS12):
            avec = np.array([[math.sin(beta) * math.sin(gamma) * apar, math.sin(
                beta) * math.cos(gamma) * apar, math.cos(beta) * apar], [
                0, bpar, -cpar], [0, bpar, cpar]], float)

        elif lattype[0:3] == b'CXZ' and abs(
                beta - math.pi / 2) < common.EPS12 and (abs(
                alpha - math.pi / 2) < common.EPS12 or abs(
                gamma - math.pi / 2) < common.EPS12):
            avec = np.array([[math.sin(gamma) * apar / 2, math.cos(
                gamma) * apar / 2, -cpar / 2], [0, math.sin(
                alpha) * bpar, math.cos(alpha) * bpar], [math.sin(
                gamma) * apar / 2, math.cos(
                gamma) * apar / 2, cpar / 2]], float)

        else:
            raise ValueError(lattype)

        return avec / apar

    def _read_file_lapw_out(self, level, filenames, lapwkunit):
        if level > 0:
            flag1, flag2 = True, True
        if level > 1:
            klist = []
        if level > 2:
            nspin = len(filenames)
            nband = sys.maxsize
            elist = []

        for ii, filename in enumerate(filenames):
            if level > 2:
                elist.append([])
            with open(filename, 'rb') as ff:
                for jj, line in enumerate(ff):
                    if ii == 0:
                        if level > 0:
                            if jj == 1:
                                self.prefix = line.strip().decode()
                            if flag1 and b'TYPE LATTICE ASSUMED' in line:
                                flag1 = False
                                lattype = line.split()[1]
                            elif flag2 and b'LATTICE CONSTANTS ARE:' in line:
                                flag2 = False
                                words = line.split()
                                apar = float(words[3])
                                bpar = float(words[4])
                                cpar = float(words[5])

                                if len(words) == 6:
                                    alpha = math.pi / 2
                                    beta = math.pi / 2
                                    gamma = math.pi / 2
                                else:
                                    alpha = float(words[6]) * math.pi / 180
                                    beta = float(words[7]) * math.pi / 180
                                    gamma = float(words[8]) * math.pi / 180

                                self.avec = self._calc_lapw_avec(
                                    lattype, apar, bpar, cpar, alpha, beta,
                                    gamma)
                                self.alat = apar
                                self.aunit = 'bohr'

                        if level > 1:
                            if b'     K=' in line:
                                dummy = np.fromstring(line[
                                    7:], dtype = float, count = 3, sep = ' ')
                                klist.append(dummy)

                    if level > 2:
                        if b'EIGENVALUES ARE:' in line:
                            buf = io.BytesIO()
                            for line in ff:
                                if b'EIGENVALUES BELOW THE ENERGY' in line:
                                    break
                                buf.write(line.strip() + b' ')
                            buf.seek(0)
                            dummy = np.loadtxt(buf, dtype = float)
                            buf.close()
                            elist[ii].append(dummy)
                            if dummy.shape[0] < nband:
                                nband = dummy.shape[0]

        if level > 1:
            self.kpoint = np.array(klist, float)
            self.kunit = lapwkunit

        if level > 2:
            nkpoint = self.kpoint.shape[0]
            energy = np.empty((nkpoint, nband, nspin), float)
            for ii in range(nspin):
                for jj in range(nkpoint):
                    energy[jj, :, ii] = elist[ii][jj][0:nband]

            self.energy = energy
            self.eunit = 'rydberg'

    def _read_file_bands_out(self, level, filenames):
        if level > 2:
            nspin = len(filenames)

        for ii, filename in enumerate(filenames):
            with open(filename, 'rb') as ff:
                line = ff.readline()
                if ii == 0:
                    words = line.split()
                    if level > 1:
                        nkpoint = int(words[4])
                        kpoint = np.empty((nkpoint, 3), float)
                    if level > 2:
                        nband = int(words[2][:-1])
                        energy = np.empty((nkpoint, nband, nspin), float)

                for jj in range(nkpoint):
                    if ii == 0 and level > 1:
                        kpoint[jj, :] = np.fromfile(
                            ff, dtype = float, count = 3, sep = ' ')
                    dummy = np.fromfile(
                        ff, dtype = float, count = nband, sep = ' ')
                    if level > 2:
                        energy[jj, :, ii] = dummy

        if level > 1:
            self.kpoint = kpoint
            self.kunit = 'cartesian'

        if level > 2:
            self.energy = energy
            self.eunit = 'ev'

    def _read_file_matdyn_out(self, level, filenames):
        if level > 1:
            klist = []
            ii = -1
        if level > 2:
            elist = []

        with open(filenames[0], 'rb') as ff:
            for line in ff:
                if level > 1:
                    if b' q =' in line:
                        dummy = np.fromstring(line[
                            4:], dtype = float, count = 3, sep = ' ')
                        klist.append(dummy)
                        ii += 1
                        elist.append([])

                if level > 2:
                    if any(ss in line for ss in [b'freq', b'omega']):
                        words = line.split()
                        elist[ii].append(float(words[7]))

        if level > 1:
            self.kpoint = np.array(klist, float)
            self.kunit = 'cartesian'

        if level > 2:
            nkpoint = len(klist)
            nband = len(elist[0])
            nspin = 1
            self.energy = np.array(elist, float).reshape(nkpoint, nband, nspin)
            self.eunit = 'cm-1'

    def _read_file_inteqp_out(self, level, filenames, etype):
        if etype is None:
            etype = 'eqp'

        if level > 2:
            dd = {'emf': 5, 'eqp': 6}
            eindex = dd[etype]

        if level > 1:
            data = np.loadtxt(filenames[0], dtype = float, skiprows = 2)
            nspin = int(round(data[-1, 0]))
            nband_exclude = int(round(data[0, 1])) - 1
            nband = int(round(data[-1, 1]))
            nkpoint = int(round(data.shape[0] / nspin / (
                    nband - nband_exclude)))

            oldaunit = self.aunit
            self.set_aunit('bohr')
            kscale = self.alat / (2 * math.pi)
            self.set_aunit(oldaunit)

            kpoint = data[0:nkpoint, 2:5] * kscale
            self.kpoint = kpoint
            self.kunit = 'cartesian'

        if level > 2:
            if nband_exclude == 0:
                energy = data[:, eindex].reshape(
                    nspin, nband, nkpoint).swapaxes(0, 2)
            else:
                energy = np.empty((nkpoint, nband, nspin), float)
                energy[:, 0:nband_exclude, :] = -common.INF12
                kk = 0
                for ii in range(nspin):
                    for jj in range(nband_exclude, nband):
                        energy[:, jj, ii] = data[kk:kk + nkpoint, eindex]
                        kk += nkpoint

            self.energy = energy
            self.eunit = 'ev'

    def _read_file_sigma_out(self, level, filenames, etype):
        if etype is None:
            etype = 'eqp1'

        if level > 1:
            klist = []
            nspin = 0
        if level > 2:
            dd = {
                'edft': 1, 'ecor': 2, 'eqp0': 8, 'eqp1': 9, 'eqp0p': 12,
                'eqp1p': 13}
            eindex = dd[etype]
            elist = []

        with open(filenames[0], 'rb') as ff:
            for line in ff:
                if level > 2:
                    if b'band_index' in line:
                        words = line.split()
                        nband = int(words[2])
                        nband_exclude = int(words[1]) - 1

                if level > 1:
                    if b'       k =' in line:
                        dummy = np.fromstring(line[
                            10:], dtype = float, count = 3, sep = ' ')
                        klist.append(dummy)

                        if level > 2:
                            ispin = int(line[line.find('spin =') + 6:])
                            if ispin > nspin:
                                nspin = ispin

                            dummy = np.genfromtxt(
                                ff, dtype = float, skip_header = 2,
                                usecols = eindex, max_rows = nband -
                                nband_exclude)
                            elist.append(dummy)

        if level > 1:
            self.kpoint = np.array(klist, float)
            self.kunit = 'crystal'

        if level > 2:
            if nband_exclude == 0:
                energy = np.array(elist, float).reshape(
                    nkpoint, nspin, nband).swapaxes(1, 2)
            else:
                energy = np.empty((nkpoint, nband, nspin), float)
                energy[:, 0:nband_exclude, :] = -common.INF12
                kk = 0
                for ii in range(nkpoint):
                    for jj in range(nspin):
                        energy[ii, nband_exclude:nband, jj] = elist[kk]
                        kk += 1

            self.energy = energy
            self.eunit = 'ev'

    def _read_file_wannier_out(self, level, filenames):
        nspin = len(filenames)

        for ii, filename in enumerate(filenames):
            if level > 1:
                data = np.loadtxt(filename, dtype = float, unpack = True)
                if ii == 0:
                    nkpoint = np.where(data[0, :] == data[0, 0])[0][1]
                    nband = data.shape[1] // nkpoint

                    kline = data[0, 0:nkpoint]
                    oldaunit = self.aunit
                    self.set_aunit('angstrom')
                    kline *= self.alat / (2 * math.pi)
                    self.kline = kline
                    self.set_aunit(oldaunit)

            if level > 2:
                if ii == 0:
                    energy = np.empty((nkpoint, nband, nspin), float)
                energy[:, :, ii] = data[1, :].reshape(nband, nkpoint).transpose(
                    )
                if ii == nspin - 1:
                    self.eunit, self.energy = 'ev', energy

    def _write_file_internal(self, level, filenames):
        s0 = '{0[0]:f} {0[1]:f} {0[2]:f}\n'
        s1 = ('{0[0][0]:2d} {0[0][1]:2d} {0[0][2]:2d}  {0[1][0]:2d} {0[1][1]:'
            '2d} {0[1][2]:2d}  {0[2][0]:2d} {0[2][1]:2d} {0[2][2]:2d}\n')
        s2 = '{0:f}\n'
        s3 = 'energy {0:d} {1:d} {2:d}\n'

        with open(filenames[0], 'wb') as ff:
            if level > 0:
                if hasattr(self, 'prefix'):
                    ff.write('prefix = {0:s}\n'.format(self.prefix).encode())
                if hasattr(self, 'aunit'):
                    ff.write('aunit = {0:s}\n'.format(self.aunit).encode())
                if hasattr(self, 'alat'):
                    ff.write('alat = {0:f}\n'.format(self.alat).encode())

                if hasattr(self, 'avec'):
                    ff.write(b'avec 3\n')
                    np.savetxt(ff, self.avec)
                if hasattr(self, 'bvec'):
                    ff.write(b'bvec 3\n')
                    np.savetxt(ff, self.bvec)

                if hasattr(self, 'avol'):
                    ff.write('avol = {0:f}\n'.format(self.avol).encode())
                if hasattr(self, 'bvol'):
                    ff.write('bvol = {0:f}\n'.format(self.bvol).encode())
                if hasattr(self, 'natom'):
                    ff.write('natom = {0:d}\n'.format(self.natom).encode())
                if hasattr(self, 'nelec'):
                    ff.write('nelec = {0:f}\n'.format(self.nelec).encode())

                if hasattr(self, 'rot'):
                    ff.write('sym {0:d}\n'.format(self.nsym).encode())
                    np.savetxt(ff, self.rot.reshape(self.nsym, 9), fmt='%d')

            if level > 1:
                if hasattr(self, 'kunit'):
                    ff.write('kunit = {0:s}\n'.format(self.kunit).encode())

                if hasattr(self, 'kpoint'):
                    ff.write('kpoint {0:d}\n'.format(self.nkpoint).encode())
                    np.savetxt(ff, self.kpoint)

                if hasattr(self, 'kline'):
                    ff.write('kline {0:d}\n'.format(self.nkpoint).encode())
                    self.kline.tofile(ff, '\n')
                    ff.write(b'\n')

                if hasattr(self, 'kweight'):
                    ff.write('kweight {0:d}\n'.format(self.nkpoint).encode())
                    self.kweight.tofile(ff, '\n')
                    ff.write(b'\n')

                if hasattr(self, 'kpath'):
                    ff.write('kpath {0:d}\n'.format(self.nkpath).encode())
                    np.savetxt(ff, self.kpath)

                if hasattr(self, 'kindex'):
                    ff.write('kindex {0:d}\n'.format(self.nkpath).encode())
                    self.kindex.tofile(ff, '\n')
                    ff.write(b'\n')

                if hasattr(self, 'klabel'):
                    ff.write('klabel {0:d}\n'.format(self.nkpath).encode())
                    for ii in range(self.nkpath):
                        ff.write('{0:s}\n'.format(self.klabel[ii]).encode())

            if level > 2:
                if hasattr(self, 'eunit'):
                    ff.write('eunit = {0:s}\n'.format(self.eunit).encode())

                if hasattr(self, 'energy'):
                    ff.write(s3.format(
                        self.nkpoint, self.nband, self.nspin).encode())
                    self.energy.tofile(ff, '\n')
                    ff.write(b'\n')

                if hasattr(self, 'efermi'):
                    ff.write('efermi = {0:f}\n'.format(self.efermi).encode())
                if hasattr(self, 'vref'):
                    ff.write('vref = {0:f}\n'.format(self.vref).encode())

    def _write_file_pw_in(self, level, filenames):
        d0 = {'bohr': 'bohr', 'angstrom': 'angstrom', 'nm': 'angstrom'}
        d1 = {'bohr': 'bohr', 'angstrom': 'angstrom'}
        d2 = {'cartesian': 'tpiba', 'crystal': 'crystal'}
        d3 = {'cartesian': 'tpiba_b', 'crystal': 'crystal_b'}

        with open(filenames[0], 'wb') as ff:
            if level > 0:
                if hasattr(self, 'aunit') and hasattr(self, 'alat') and hasattr(
                        self, 'avec'):
                    oldaunit = self.aunit
                    self.set_aunit(d0[self.aunit])

                    ff.write('CELL_PARAMETERS {0:s}\n'.format(d1[
                        self.aunit]).encode())
                    np.savetxt(ff, self.alat * self.avec)
                    self.set_aunit(oldaunit)

            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpoint'):
                    ff.write('K_POINTS {0:s}\n'.format(d2[self.kunit]).encode())
                    ff.write('{0:d}\n'.format(self.nkpoint).encode())

                    if hasattr(self, 'kweight'):
                        weight = 2.0 * self.kweight
                    else:
                        weight = np.full(self.nkpoint, 2 / self.nkpoint)
                    np.savetxt(ff, np.concatenate((
                        self.kpoint, weight.reshape(self.nkpoint, 1)), 1))

                if hasattr(self, 'kunit') and hasattr(
                        self, 'kpath') and hasattr(self, 'kindex'):
                    if self.check_kindex() == 'uninitialized':
                        raise ValueError('call calc_kindex')

                    ff.write('K_POINTS {0:s}\n'.format(d3[self.kunit]).encode())
                    ff.write('{0:d}\n'.format(self.nkpath).encode())

                    knum = np.append(np.diff(self.kindex), 0)
                    np.savetxt(ff, np.concatenate((
                        self.kpath, knum.reshape(self.nkpath, 1)), 1))

    def _write_file_wannier_in(self, level, filenames):
        d0 = {'bohr': 'Bohr', 'angstrom': 'Angstrom'}
        l0 = [' ', '\n']

        with open(filenames[0], 'wb') as ff:
            if level > 0:
                if hasattr(self, 'aunit') and hasattr(self, 'alat') and hasattr(
                        self, 'avec'):
                    ff.write(b'Begin Unit_Cell_Cart\n')
                    ff.write('{0:s}'.format(d0[self.aunit]).encode())
                    np.savetxt(ff, self.alat * self.avec)
                    ff.write(b'End Unit_Cell_Cart\n')

            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                    oldkunit = self.kunit
                    self.set_kunit('crystal')

                    ff.write(b'begin kpoint_path\n')
                    for ii in range(self.nkpath - 1):
                        for jj in range(2):
                            ff.write('{0:s} '.format(self.klabel[
                                ii + jj]).encode())
                            self.kpath[ii + jj, :].tofile(ff, ' ')
                            ff.write('{0:s}'.format(l0[jj]).encode())
                    ff.write(b'end kpoint_path\n')

                    self.set_kunit(oldkunit)

    def _write_file_vasp_kpt(self, level, filenames):
        d0 = {'cartesian': 'cartesian', 'crystal': 'reciprocal'}

        with open(filenames[0], 'wb') as ff:
            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                    if self.check_kindex() != 'number':
                        raise ValueError("call calc_kindex('number')")

                    ff.write(b'k-points along high symmetry lines\n')
                    ff.write('{0:d}\n'.format(self.kindex[1]).encode())
                    ff.write(b'line-mode\n')
                    ff.write('{0:s}'.format(d0[self.kunit]).encode())

                    for ii in range(self.nkpath - 1):
                        ff.write(b'\n')
                        for jj in range(2):
                            self.kpath[ii + jj, :].tofile(ff, ' ')
                            ff.write(' ! {0:s}\n'.format(self.klabel[
                                ii + jj]).encode())

    def _write_file_lapw_kpt(self, level, filenames, lapwkunit):
        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                if self.check_kindex() == 'uninitialized':
                    raise ValueError('call calc_kindex')

                oldkunit = self.kunit
                self.set_kunit(lapwkunit)
                nline = self.nkpath - 1
                nstep = np.diff(self.kindex)
                line = np.diff(self.kpath, axis = 0)
                step = line / np.broadcast_to(nstep, (3, nline)).transpose()
                origin = self.kpath[0, :] - step[0, :]
                self.set_kunit(oldkunit)

                ngridmax = 8192
                ngrid = []
                for ii in range(nline):
                    denom = []
                    for jj in range(3):
                        aa = fractions.Fraction(step[ii, jj]).limit_denominator(
                            )
                        if aa.numerator != 0 and aa.denominator not in denom:
                            denom.append(aa.denominator)

                    ll = denom[0]
                    for jj in range(1, len(denom)):
                        ll *= denom[jj] // fractions.gcd(ll, denom[jj])
                    ngrid.append(min(ll, ngridmax))

                labels = []
                kpoints = []
                ngrids = []
                weights = []
                kvec = origin
                for ii in range(nline):
                    for jj in range(min(ii, 1), nstep[ii] + 1):
                        if jj == 0:
                            label = '{0:<10s}'.format(self.klabel[0])
                        elif jj == nstep[ii]:
                            label = '{0:<10s}'.format(self.klabel[ii + 1])
                        else:
                            label = ' ' * 10

                        weight = '  2.0'
                        if jj == 0:
                            weight += ' -1.0 1.5'
                            if hasattr(self, 'prefix'):
                                weight += '      {0:s}'.format(self.prefix)

                        kvec += step[ii]
                        ix = int(round(kvec[0] * ngrid[ii]))
                        iy = int(round(kvec[1] * ngrid[ii]))
                        iz = int(round(kvec[2] * ngrid[ii]))

                        labels.append(label)
                        kpoints.append([ix, iy, iz])
                        ngrids.append(ngrid[ii])
                        weights.append(weight)

                s0 = '{0:s}{1[0]:5d}{1[1]:5d}{1[2]:5d}{2:5d}{3:s}\n'

                with open(filenames[0], 'wb') as ff:
                    for ii in range(len(labels)):
                        ff.write(s0.format(labels[ii], kpoints[ii], ngrids[
                            ii], weights[ii]).encode())
                    ff.write(b'END\n')

    def _write_file_boltztrap_in(self, level, filenames, boltzparam):
        if level > 2:
            self.set_aunit('bohr')
            self.set_alat(1.0)
            self.set_kunit('crystal')
            self.set_eunit('rydberg')

            (deltae, ecut, lpfac, efcut, tmax, deltat, ecut2, dosmethod,
                nband_exclude) = boltzparam

            with open(filenames[0], 'wb') as ff:
                filename_intrans = os.path.basename(filenames[1])
                filename_struct = os.path.basename(filenames[2])
                filename_energy = os.path.basename(filenames[3])

                s0 = ", 'old', 'formatted', 0\n"
                s1 = ", 'unknown', 'formatted', 0\n"

                ff.write("5, '{0:s}'{1:s}".format(
                    filename_intrans, s0).encode())
                ff.write("6, '{0:s}.outputtrans'{1:s}".format(
                    self.prefix, s1).encode())
                ff.write("20, '{0:s}'{1:s}".format(
                    filename_struct, s0).encode())
                ff.write("10, '{0:s}'{1:s}".format(
                    filename_energy, s0).encode())
                ff.write("88, '{0:s}.epa.e'{1:s}".format(
                    self.prefix, s0).encode())
                ff.write("-1, '{0:s}'{1:s}".format(
                    self.prefix, s0).encode())

            with open(filenames[1], 'wb') as ff:
                nelec = self.nelec - 2 * nband_exclude

                s0 = '# Format of DOS\n'
                s1 = '# iskip (not presently used) idebug setgap shiftgap\n'
                s2 = ('# Fermilevel (Ry), energygrid, energy span around '
                    'Fermilevel, number of electrons\n')
                s3 = ('# CALC (calculate expansion coeff), NOCALC read from '
                    'file\n')
                s4 = '# lpfac, number of latt-points per k-point\n'
                s5 = '# run mode (only BOLTZ is supported)\n'
                s6 = '# (efcut) energy range of chemical potential\n'
                s7 = '# Tmax, temperature grid\n'
                s8 = ('# energyrange of bands given individual DOS output '
                    'sig_xxx and dos_xxx (xxx is band number)\n')

                ff.write('GENE      {0:s}'.format(s0).encode())
                ff.write('0 0 0 0.0 {0:s}'.format(s1).encode())
                ff.write('{0:f} {1:f} {2:f} {3:f} {4:s}'.format(
                    self.efermi, deltae, ecut, nelec, s2).encode())
                ff.write('CALC {0:s}'.format(s3).encode())
                ff.write('{0:d} {1:s}'.format(lpfac, s4).encode())
                ff.write('BOLTZ {0:s}'.format(s5).encode())
                ff.write('{0:f} {1:s}'.format(efcut, s6).encode())
                ff.write('{0:f} {1:f} {2:s}'.format(tmax, deltat, s7).encode())
                ff.write('{0:f} {1:s}'.format(ecut2, s8).encode())
                ff.write('{0:s}\n'.format(dosmethod).encode())

            with open(filenames[2], 'wb') as ff:
                ff.write('{0:s}\n'.format(self.prefix).encode())
                np.savetxt(ff, self.avec)
                ff.write('{0:d}\n'.format(self.nsym).encode())
                np.savetxt(ff, self.rot.transpose(0, 2, 1).reshape(
                    self.nsym, 9), fmt='%d')

            with open(filenames[3], 'wb') as ff:
                nband = self.nspin * (self.nband - nband_exclude)
                ff.write('{0:s}\n'.format(self.prefix).encode())
                ff.write('{0:d}\n'.format(self.nkpoint).encode())
                for ii in range(self.nkpoint):
                    self.kpoint[ii, :].tofile(ff, ' ')
                    ff.write(' {0:d}\n'.format(nband).encode())
                    self.energy[ii, nband_exclude:self.nband, :].tofile(
                        ff, '\n')
                    ff.write(b'\n')

