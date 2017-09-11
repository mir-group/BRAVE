"""This module defines class File."""

import sys
import math
import fractions
import os

import numpy

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
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 0:
            iprefix = -1
            iaunit = -1
            ialat = -1
            iavec = -1
            ibvec = -1
            iavol = -1
            ibvol = -1
            inatom = -1
            inelec = -1
            isym = -1
        if level > 1:
            ikunit = -1
            ikpt = -1
            ikline = -1
            ikweight = -1
            ikpath = -1
            ikindex = -1
            iklabel = -1
        if level > 2:
            ieunit = -1
            ienergy = -1
            iefrm = -1
            ivref = -1
        for ii, line in enumerate(contents[0]):
            if line.lstrip()[0:1] == '#':
                continue
            if level > 0:
                if 'prefix' in line:
                    iprefix = ii
                elif 'aunit' in line:
                    iaunit = ii
                elif 'alat' in line:
                    ialat = ii
                elif 'avec' in line:
                    iavec = ii
                elif 'bvec' in line:
                    ibvec = ii
                elif 'avol' in line:
                    iavol = ii
                elif 'bvol' in line:
                    ibvol = ii
                elif 'natom' in line:
                    inatom = ii
                elif 'nelec' in line:
                    inelec = ii
                elif 'sym' in line:
                    isym = ii
            if level > 1:
                if 'kunit' in line:
                    ikunit = ii
                elif 'kpoint' in line:
                    ikpt = ii
                elif 'kline' in line:
                    ikline = ii
                elif 'kweight' in line:
                    ikweight = ii
                elif 'kpath' in line:
                    ikpath = ii
                elif 'kindex' in line:
                    ikindex = ii
                elif 'klabel' in line:
                    iklabel = ii
            if level > 2:
                if 'eunit' in line:
                    ieunit = ii
                elif 'energy' in line:
                    ienergy = ii
                elif 'efermi' in line:
                    iefrm = ii
                elif 'vref' in line:
                    ivref = ii

        if level > 0:
            if iprefix != -1:
                self.prefix = contents[0][iprefix].replace('=', ' = ').split(
                        )[2]
            if iaunit != -1:
                self.aunit = contents[0][iaunit].replace('=', ' = ').split()[2]
            if ialat != -1:
                self.alat = float(contents[0][ialat].replace('=', ' = ').split(
                        )[2])

            if iavec != -1:
                adim = int(contents[0][iavec].split()[1])
                if adim != 3:
                    raise ValueError(adim)

                avec = []
                for ii in range(3):
                    words = contents[0][iavec + 1 + ii].split()
                    avec.append([])
                    for word in words:
                        avec[ii].append(float(word))
                self.avec = avec

            if ibvec != -1:
                bdim = int(contents[0][ibvec].split()[1])
                if bdim != 3:
                    raise ValueError(bdim)

                bvec = []
                for ii in range(3):
                    words = contents[0][ibvec + 1 + ii].split()
                    bvec.append([])
                    for word in words:
                        bvec[ii].append(float(word))
                self.bvec = bvec

            if iavol != -1:
                self.avol = float(contents[0][iavol].replace('=', ' = ').split(
                        )[2])
            if ibvol != -1:
                self.bvol = float(contents[0][ibvol].replace('=', ' = ').split(
                        )[2])
            if inatom != -1:
                self.natom = int(contents[0][inatom].replace('=', ' = ').split(
                        )[2])
            if inelec != -1:
                self.nelec = float(contents[0][inelec].replace(
                        '=', ' = ').split()[2])

            if isym != -1:
                nsym = int(contents[0][isym].split()[1])
                rot = numpy.empty((nsym, 3, 3), int)
                for ii in range(nsym):
                    words = contents[0][isym + 1 + ii].split()
                    for jj in range(3):
                        for kk in range(3):
                            rot[ii, jj, kk] = float(words[jj * 3 + kk])
                self.rot = rot

            if iavec != -1 and ibvec == -1:
                self.calc_bvec()
            elif iavec == -1 and ibvec != -1:
                self.calc_avec()

        if level > 1:
            if ikunit != -1:
                self.kunit = contents[0][ikunit].replace('=', ' = ').split()[2]

            if ikpt != -1:
                nkpoint = int(contents[0][ikpt].split()[1])
                kpoint = numpy.empty((nkpoint, 3), float)
                for ikpoint in range(nkpoint):
                    words = contents[0][ikpt + 1 + ikpoint].split()
                    for jj in range(3):
                        kpoint[ikpoint, jj] = float(words[jj])
                self.kpoint = kpoint

            if ikline != -1:
                nkpoint = int(contents[0][ikline].split()[1])
                kline = numpy.empty(nkpoint, float)
                for ikpoint in range(nkpoint):
                    kline[ikpoint] = float(contents[0][ikline + 1 + ikpoint])
                self.kline = kline

            if ikweight != -1:
                nkpoint = int(contents[0][ikweight].split()[1])
                kweight = numpy.empty(nkpoint, float)
                for ikpoint in range(nkpoint):
                    kweight[ikpoint] = float(contents[0][
                            ikweight + 1 + ikpoint])
                self.kweight = kweight

            if ikpath != -1:
                nkpath = int(contents[0][ikpath].split()[1])
                kpath = numpy.empty((nkpath, 3), float)
                for ikpoint in range(nkpath):
                    words = contents[0][ikpath + 1 + ikpoint].split()
                    for jj in range(3):
                        kpath[ikpoint, jj] = float(words[jj])
                self.kpath = kpath

            if ikindex != -1:
                nkpath = int(contents[0][ikindex].split()[1])
                kindex = numpy.empty(nkpath, int)
                for ikpoint in range(nkpath):
                    kindex[ikpoint] = int(contents[0][ikindex + 1 + ikpoint])
                self.kindex = kindex

            if iklabel != -1:
                nkpath = int(contents[0][iklabel].split()[1])
                klabel = []
                for ikpoint in range(nkpath):
                    words = contents[0][iklabel + 1 + ikpoint].split()
                    klabel.append(words[0])
                self.klabel = klabel

        if level > 2:
            if ieunit != -1:
                self.eunit = contents[0][ieunit].replace('=', ' = ').split()[2]

            if ienergy != -1:
                nkpoint = int(contents[0][ienergy].split()[1])
                if nkpoint != self.nkpoint:
                    raise ValueError(nkpoint)
                nband = int(contents[0][ienergy].split()[2])
                nspin = int(contents[0][ienergy].split()[3])
                energy = numpy.empty((nkpoint, nband, nspin), float)
                ii = ienergy + 1
                for ikpoint in range(nkpoint):
                    for iband in range(nband):
                        for ispin in range(nspin):
                            energy[ikpoint, iband, ispin] = float(
                                    contents[0][ii].split()[0])
                            ii += 1
                self.energy = energy

            if iefrm != -1:
                self.efermi = float(contents[0][iefrm].replace(
                        '=', ' = ').split()[2])

            if ivref != -1:
                self.vref = float(contents[0][ivref].replace('=', ' = ').split(
                        )[2])

    def _read_file_pw_out(self, level, filenames):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 0:
            avec = []
            bvec = []
            idxsym = []
        if level > 1:
            nspin = 1
        if level > 2:
            idxbnd = []
            iefermi = 0
        for ii, line in enumerate(contents[0]):
            if level > 0:
                if 'Writing output data file' in line:
                    prefix = line.split()[4][:-5]

                elif 'lattice parameter (alat)  =' in line or (
                        'lattice parameter (a_0)   =' in line):
                    alat = float(line.split()[4])

                elif ('a(1) = (' in line or 'a(2) = (' in line
                        or 'a(3) = (' in line):
                    words = line.replace('(', ' ').replace(')', ' ').split()
                    avec.append((
                            float(words[3]), float(words[4]), float(words[5])))

                elif ('b(1) = (' in line or 'b(2) = (' in line
                        or 'b(3) = (' in line):
                    words = line.replace('(', ' ').replace(')', ' ').split()
                    bvec.append((
                            float(words[3]), float(words[4]), float(words[5])))

                elif 'unit-cell volume          =' in line:
                    avol = float(line.split()[3])

                elif 'number of atoms/cell' in line:
                    natom = int(line.split()[4])

                elif 'number of electrons' in line:
                    nelec = float(line.split()[4])

                elif 'cryst.   s' in line:
                    idxsym.append(ii)

            if level > 1:
                if 'number of k points=' in line:
                    nkpoint = int(line.split()[4])
                elif 'SPIN' in line:
                    nspin = 2
                elif ' cryst. coord.' in line:
                    idxkpt = ii + 1

            if level > 2:
                if 'number of Kohn-Sham states=' in line:
                    nband = int(line.split()[4])
                elif 'band energies (ev)' in line or 'bands (ev)' in line:
                    idxbnd.append(ii + 2)
                elif 'the Fermi energy is' in line:
                    efermi = float(line.split()[4])
                    iefermi += 1
                elif 'highest occupied, lowest unoccupied level' in line:
                    efermi = (float(line.split()[6]) + float(
                            line.split()[7])) / 2.0
                    iefermi += 1

        if level > 0:
            nsym = len(idxsym)
            if nsym > 0:
                rot = numpy.empty((nsym, 3, 3), int)
                for ii in range(nsym):
                    for jj in range(3):
                        words = contents[0][idxsym[ii] + jj][19:53].split()
                        for kk in range(3):
                            rot[ii, jj, kk] = int(words[kk])

            self.prefix, self.natom, self.nelec = prefix, natom, nelec
            self.aunit, self.alat, self.avec, self.bvec, self.avol = (
                    'bohr', alat, avec, bvec, avol / alat ** 3)
            if nsym > 0:
                self.rot = rot

        if level > 1:
            nkpoint //= nspin
            kpoint = numpy.empty((nkpoint, 3), float)
            kweight = numpy.empty(nkpoint, float)
            for ikpoint in range(nkpoint):
                words = contents[0][idxkpt + ikpoint][20:56].split()
                words.append(contents[0][idxkpt + ikpoint].split()[-1])
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj])
                kweight[ikpoint] = 0.5 * float(words[3])

            self.kunit, self.kpoint, self.kweight = 'crystal', kpoint, kweight

        if level > 2:
            energy = numpy.empty((nkpoint, nband, nspin), float)

            ncol = 8
            nrow = nband // ncol
            if nband % ncol != 0:
                nrow += 1

            for ikpoint in range(nkpoint):
                for ispin in range(nspin):
                    iband = 0
                    for irow in range(nrow):
# This does not work if there are no spaces between eigenvalues
# which often happens for semicore states.
#                        words = contents[0][idxbnd[
#                                ikpoint + ispin * nkpoint] + irow].split()
# This works even if there are no spaces between eigenvalues.
                        line = contents[0][idxbnd[
                                ikpoint + ispin * nkpoint] + irow]
                        words = []
                        for jj in range(2, len(line) - 9, 9):
                            words.append(line[jj:jj + 9])
                        for word in words:
                            energy[ikpoint, iband, ispin] = float(word)
                            iband += 1

            self.eunit, self.energy = 'ev', energy
            if iefermi != 0:
                self.efermi = efermi

    def _read_file_wannier_in(self, level, filenames):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 1:
            nkpersect = 100
        for ii, line in enumerate(contents[0]):
            lline = line.lower()
            if level > 0:
                if 'begin' in lline and 'unit_cell_cart' in lline:
                    iaunit = ii + 1
            if level > 1:
                if 'bands_num_points' in lline:
                    nkpersect = int(line.split()[2])
                elif 'begin' in lline and 'kpoint_path' in lline:
                    istart = ii
                elif 'end' in lline and 'kpoint_path' in lline:
                    iend = ii

        if level > 0:
            if contents[0][iaunit].split()[0].lower() == 'bohr':
                aunit = 'bohr'
            else:
                aunit = 'angstrom'

            avec = []
            for ii in range(3):
                words = contents[0][iaunit + 1 + ii].split()
                avec.append([])
                for word in words:
                    avec[ii].append(float(word))

            self.aunit, self.alat, self.avec = aunit, 1, avec
            self.calc_bvec()

        if level > 1:
            nkpath = iend - istart

            k1label = []
            k1coord = []
            k2label = []
            k2coord = []
            for ii in range(nkpath - 1):
                words = contents[0][istart + 1 + ii].split()
                k1label.append(words[0])
                k1coord.append([])
                k2label.append(words[4])
                k2coord.append([])
                for jj in range(3):
                    k1coord[ii].append(float(words[1 + jj]))
                    k2coord[ii].append(float(words[5 + jj]))

            if k1label[1:] != k2label[:-1] or k1coord[1:] != k2coord[:-1]:
                raise ValueError(kpath)

            kpath = k1coord + [k2coord[nkpath - 1]]
            kindex = [nkpersect] + [0 for ii in range(nkpath - 1)]
            klabel = k1label + [k2label[nkpath - 1]]

            self.kunit, self.kpath, self.kindex, self.klabel = (
                    'crystal', kpath, kindex, klabel)

    def _read_file_vasp_out(self, level, filenames):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 0:
            avec = []
        if level > 2:
            nspin = 1
            idxbnd = []
            iefermi = 0
        for ii, line in enumerate(contents[0]):
            if level > 0:
                if 'SYSTEM =' in line:
                    prefix = line.split()[2]

                if 'ALAT       =' in line:
                    alat = float(line.split()[2])

                elif 'A1 = (' in line or 'A2 = (' in line or 'A3 = (' in line:
                    words = line.replace(',', ' ').replace('(', ' ').replace(
                            ')', ' ').split()
                    avec.append((
                            float(words[2]), float(words[3]), float(words[4])))
            if level > 1:
                if 'k-points           NKPTS =' in line:
                    nkpoint = int(line.split()[3])

                elif 'k-points in reciprocal lattice and weights:' in line:
                    idxkpt = ii + 1

            if level > 2:
                if 'k-points NKPTS =' in line and (
                        'number of bands NBANDS=' in line):
                    nkpoint = int(line.split()[3])
                    if nkpoint != self.nkpoint:
                        raise ValueError(nkpoint)
                    nband = int(line.split()[14])

                elif 'ISPIN' in line:
                    nspin = int(line.split()[2])

                elif 'band No.  band energies occupation' in line:
                    idxbnd.append(ii + 1)

                elif 'E-fermi :' in line:
                    efermi = float(line.split()[2])
                    iefermi += 1

        if level > 0:
            self.prefix, self.aunit, self.alat, self.avec = (
                    prefix, 'angstrom', 1, avec)
            self.calc_bvec()
            self.set_alat(alat)

        if level > 1:
            kpoint = numpy.empty((nkpoint, 3), float)
            for ikpoint in range(nkpoint):
                words = contents[0][idxkpt + ikpoint].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj])

            self.kunit, self.kpoint = 'crystal', kpoint

        if level > 2:
            energy = numpy.empty((nkpoint, nband, nspin), float)
            ii = 0
            for ispin in range(nspin):
                for ikpoint in range(nkpoint):
                    for iband in range(nband):
                        words = contents[0][idxbnd[ii] + iband].split()
                        energy[ikpoint, iband, ispin] = float(words[1])
                    ii += 1

            self.eunit, self.energy = 'ev', energy
            if iefermi != 0:
                self.efermi = efermi

    def _read_file_lapw_out(self, level, filenames, lapwkunit):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 0:
            avec = []
            prefix = contents[0][1].split()[0]
        if level > 1:
            idxkpt = []
        if level > 2:
            nspin = len(contents)
            idxbnd1 = []
            idxbnd2 = []
        for ispin in range(nspin):
            if level > 2:
                idxbnd1.append([])
                idxbnd2.append([])
            for ii, line in enumerate(contents[ispin]):
                if level > 0 and ispin == 0:
                    if 'TYPE LATTICE ASSUMED' in line:
                        lattice = line.split()[1]

                    elif 'LATTICE CONSTANTS ARE:' in line:
                        words = line.split()
                        apar = float(words[3])
                        bpar = float(words[4])
                        cpar = float(words[5])

                        if len(words) == 6:
                            alpha = 90.0 * math.pi / 180.0
                            beta = 90.0 * math.pi / 180.0
                            gamma = 90.0 * math.pi / 180.0
                        else:
                            alpha = float(words[6]) * math.pi / 180.0
                            beta = float(words[7]) * math.pi / 180.0
                            gamma = float(words[8]) * math.pi / 180.0

                if level > 1 and ispin == 0:
                    if '     K=' in line:
                        idxkpt.append(ii)

                if level > 2:
                    if 'EIGENVALUES ARE:' in line:
                        idxbnd1[ispin].append(ii + 1)
                    elif 'EIGENVALUES BELOW THE ENERGY' in line:
                        idxbnd2[ispin].append(ii)

        if level > 0:
            if lattice[0:1] == 'P':
                phi = math.acos((math.cos(gamma) - math.cos(alpha) * math.cos(
                        beta)) / (math.sin(alpha) * math.sin(beta)))
                avec.append((math.sin(beta) * math.sin(phi) * apar, math.sin(
                        beta) * math.cos(phi) * apar, math.cos(beta) * apar))
                avec.append((0.0, math.sin(alpha) * bpar, math.cos(
                        alpha) * bpar))
                avec.append((0.0, 0.0, cpar))

            elif lattice[0:1] == 'H':
                avec.append((0.5 * math.sqrt(3.0) * apar, -0.5 * apar, 0.0))
                avec.append((0.0, apar, 0.0))
                avec.append((0.0, 0.0, cpar))

            elif lattice[0:1] == 'T' or lattice[0:1] == 'R':
                avec.append(((0.5 / math.sqrt(3.0)) * apar, -0.5 * apar, (
                        1.0 / 3.0) * cpar))
                avec.append(((0.5 / math.sqrt(3.0)) * apar, 0.5 * apar, (
                        1.0 / 3.0) * cpar))
                avec.append(((-1.0 / math.sqrt(3.0)) * apar, 0.0, (
                        1.0 / 3.0) * cpar))

            elif lattice[0:1] == 'F':
                avec.append((0.0, 0.5 * bpar, 0.5 * cpar))
                avec.append((0.5 * apar, 0.0, 0.5 * cpar))
                avec.append((0.5 * apar, 0.5 * bpar, 0.0))

            elif lattice[0:1] == 'B':
                avec.append((-0.5 * apar, 0.5 * bpar, 0.5 * cpar))
                avec.append((0.5 * apar, -0.5 * bpar, 0.5 * cpar))
                avec.append((0.5 * apar, 0.5 * bpar, -0.5 * cpar))

            elif lattice[0:3] == 'CXY' and abs(
                    gamma - 0.5 * math.pi) < common.EPS12 and (abs(
                    alpha - 0.5 * math.pi) < common.EPS12 or abs(
                    beta - 0.5 * math.pi) < common.EPS12):
                avec.append((0.5 * math.sin(beta) * apar, -0.5 * math.sin(
                        alpha) * bpar, 0.5 * math.cos(
                        beta) * apar - 0.5 * math.cos(alpha) * bpar))
                avec.append((0.5 * math.sin(beta) * apar, 0.5 * math.sin(
                        alpha) * bpar, 0.5 * math.cos(
                        beta) * apar + 0.5 * math.cos(alpha) * bpar))
                avec.append((0.0, 0.0, cpar))

            elif lattice[0:3] == 'CYZ' and abs(
                    alpha - 0.5 * math.pi) < common.EPS12 and (abs(
                    beta - 0.5 * math.pi) < common.EPS12 or abs(
                    gamma - 0.5 * math.pi) < common.EPS12):
                avec.append((math.sin(beta) * math.sin(gamma) * apar, math.sin(
                        beta) * math.cos(gamma) * apar, math.cos(beta) * apar))
                avec.append((0.0, bpar, -cpar))
                avec.append((0.0, bpar, cpar))

            elif lattice[0:3] == 'CXZ' and abs(
                    beta - 0.5 * math.pi) < common.EPS12 and (abs(
                    alpha - 0.5 * math.pi) < common.EPS12 or abs(
                    gamma - 0.5 * math.pi) < common.EPS12):
                avec.append((0.5 * math.sin(gamma) * apar, 0.5 * math.cos(
                        gamma) * apar, -0.5 * cpar))
                avec.append((0.0, math.sin(alpha) * bpar, math.cos(
                        alpha) * bpar))
                avec.append((0.5 * math.sin(gamma) * apar, 0.5 * math.cos(
                        gamma) * apar, 0.5 * cpar))

            else:
                raise ValueError(lattice)

            self.prefix, self.aunit, self.alat, self.avec, self.kunit = (
                    prefix, 'bohr', 1, avec, lapwkunit)
            self.calc_bvec()
            self.set_alat(apar)

        if level > 1:
            nkpoint = len(idxkpt)

            kpoint = numpy.empty((nkpoint, 3), float)
            for ikpoint in range(nkpoint):
                words = contents[0][idxkpt[ikpoint]].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[1 + jj])

            self.kpoint = kpoint

        if level > 2:
            nband = sys.maxsize
            for ispin in range(nspin):
                for ikpoint in range(nkpoint):
                    iband = 0
                    for irow in range(
                            idxbnd1[ispin][ikpoint], idxbnd2[ispin][ikpoint]):
                        iband += len(contents[ispin][irow].split())
                    if iband < nband:
                        nband = iband

            energy = numpy.zeros((nkpoint, nband, nspin), float)
            for ispin in range(nspin):
                for ikpoint in range(nkpoint):
                    iband = 0
                    for irow in range(
                            idxbnd1[ispin][ikpoint], idxbnd2[ispin][ikpoint]):
                        words = contents[ispin][irow].split()
                        for word in words:
                            if iband < nband:
                                energy[ikpoint, iband, ispin] = float(word)
                            iband += 1

            self.eunit, self.energy = 'rydberg', energy

    def _read_file_bands_out(self, level, filenames):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 1:
            nkpoint = int(contents[0][0].split()[4])
            nband = int(contents[0][0].split()[2].replace(',', ''))
            ncol = 10
            nrow = nband // ncol
            if nband % ncol != 0:
                nrow += 1

            kpoint = numpy.empty((nkpoint, 3), float)
            for ikpoint in range(nkpoint):
                words = contents[0][ikpoint * (nrow + 1) + 1].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj])

            self.kunit, self.kpoint = 'cartesian', kpoint

        if level > 2:
            nspin = len(contents)

            energy = numpy.empty((nkpoint, nband, nspin), float)
            for ispin in range(nspin):
                for ikpoint in range(nkpoint):
                    iband = 0
                    for irow in range(nrow):
                        words = contents[ispin][
                                ikpoint * (nrow + 1) + irow + 2].split()
                        for word in words:
                            energy[ikpoint, iband, ispin] = float(word)
                            iband += 1

            self.eunit, self.energy = 'ev', energy

    def _read_file_matdyn_out(self, level, filenames):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 1:
            idxkpt = []
        if level > 2:
            idxbnd = []
        ii = 0
        for line in contents[0]:
            if level > 1:
                if ' q =' in line:
                    idxkpt.append(ii)
            if level > 2:
                if 'freq' in line or 'omega' in line:
                    idxbnd.append(ii)

        if level > 1:
            nkpoint = len(idxkpt)
            kpoint = numpy.empty((nkpoint, 3), float)
            for ikpoint in range(nkpoint):
                words = contents[0][idxkpt[ikpoint]].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj + 2])

            self.kunit, self.kpoint = 'cartesian', kpoint

        if level > 2:
            nband = len(idxbnd) // nkpoint
            nspin = 1

            energy = numpy.empty((nkpoint, nband, nspin), float)
            ispin = 0
            for ikpoint in range(nkpoint):
                for iband in range(nband):
                    line = contents[0][idxbnd[ikpoint*nband+iband]]
                    idx1 = len(line)-line[::-1].find('=')
                    idx2 = len(line)-line[::-1].find('[')
                    energy[ikpoint, iband, ispin] = float(line[idx1:idx2 - 1])

            self.eunit, self.energy = 'cm-1', energy

    def _read_file_inteqp_out(self, level, filenames, etype):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 1:
            oldaunit = self.aunit
            self.set_aunit('bohr')
            kscale = self.alat / (2.0 * math.pi)

        if level > 2:
            if etype is None:
                etype = 'eqp'
            dd = {'emf': 5, 'eqp': 6}
            eindex = dd[etype]

        if level > 1:
            nhead = 2
            ndata = len(contents[0]) - nhead
            nspin = int(contents[0][nhead + ndata - 1].split()[0])
            nfirst = int(contents[0][nhead].split()[1])
            nband = int(contents[0][nhead + ndata - 1].split()[1])
            nkpoint = ndata // ((nband - nfirst + 1) * nspin)

            kpoint = numpy.empty((nkpoint, 3), float)
            for ikpoint in range(nkpoint):
                words = contents[0][nhead + ikpoint][14:50].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj]) * kscale

            self.kunit, self.kpoint = 'cartesian', kpoint
            self.set_aunit(oldaunit)

        if level > 2:
            energy = numpy.empty((nkpoint, nband, nspin), float)
            energy[:, 0:nfirst - 1, :] = -common.INF12
            ii = nhead
            for ispin in range(nspin):
                for iband in range(nfirst - 1, nband):
                    for ikpoint in range(nkpoint):
                        energy[ikpoint, iband, ispin] = float(
                                contents[ii].split()[eindex])
                        ii += 1

            self.eunit, self.energy = 'ev', energy

    def _read_file_sigma_out(self, level, filenames, etype):
        contents = []
        for filename in filenames:
            with open(filename) as fileobj:
                content = fileobj.readlines()
            contents.append(content)

        if level > 2:
            if etype is None:
                etype = 'eqp1'
            dd = {
                    'edft': 1, 'ecor': 2, 'eqp0': 8, 'eqp1': 9, 'eqp0p': 12,
                    'eqp1p': 13}
            eindex = dd[etype]

        if level > 1:
            nkpoint = 0
            idxkpt = []
            nspin = 0
        if level > 2:
            nband = 0
        for ii, line in enumerate(contents[0]):
            if level > 1:
                if '       k =' in line:
                    nkpoint += 1
                    idxkpt.append(ii)
                    words = line.split()
                    nspin = max(nspin, int(words[10]))
            if level > 2:
                if 'band_index' in line:
                    words = line.split()
                    nband = int(words[2]) - int(words[1]) + 1

        if level > 1:
            kpoint = numpy.empty((nkpoint, 3), float)

            for ikpoint in range(nkpoint):
                words = contents[0][idxkpt[ikpoint]].split()
                for jj in range(3):
                    kpoint[ikpoint, jj] = float(words[jj + 2])

            self.kunit, self.kpoint = 'crystal', kpoint

        if level > 2:
            energy = numpy.empty((nkpoint, nband, nspin), float)

            for ikpoint in range(nkpoint):
                for ispin in range(nspin):
                    for iband in range(nband):
                        energy[ikpoint, iband, ispin] = float(
                                contents[0][idxkpt[ikpoint * nspin + ispin
                                ] + 3 + iband].split()[eindex])

            self.eunit, self.energy = 'ev', energy

    def _read_file_wannier_out(self, level, filenames):
        nspin = len(filenames)
        for ii, filename in enumerate(finenames):
            if level > 1:
                wan = numpy.loadtxt(filename, dtype = float, unpack = True)
                if ii == 0:
                    nkpoint = numpy.where(wan[0, :] == wan[0, 0])[0][1]
                    nband = wan.shape[0] // nkpoint

                    kline = wan[0, 0:nkpoint]
                    oldaunit = self.aunit
                    self.set_aunit('angstrom')
                    kline *= self.alat / (2.0 * math.pi)
                    self.kline = kline
                    self.set_aunit(oldaunit)

            if level > 2:
                if ii == 0:
                    energy = numpy.empty((nkpoint, nband, nspin), float)
                energy[:, :, ii] = wan[1, :].reshape(nband, nkpoint).transpose()
                self.eunit, self.energy = 'ev', energy

    def _write_file_internal(self, level, filenames):
        s0 = '{0[0]:f} {0[1]:f} {0[2]:f}\n'
        s1 = '{0[0][0]:2d} {0[0][1]:2d} {0[0][2]:2d}  {0[1][0]:2d} {0[1][1]:2d}'
        s1 += ' {0[1][2]:2d}  {0[2][0]:2d} {0[2][1]:2d} {0[2][2]:2d}\n'
        s2 = '{0:f}\n'
        s3 = 'energy {0:d} {1:d} {2:d}\n'

        with open(filenames[0], 'wb') as ff:
            if level > 0:
                if hasattr(self, 'prefix'):
                    ff.write(b'prefix = {0:s}\n'.format(self.prefix))
                if hasattr(self, 'aunit'):
                    ff.write(b'aunit = {0:s}\n'.format(self.aunit))
                if hasattr(self, 'alat'):
                    ff.write(b'alat = {0:f}\n'.format(self.alat))

                if hasattr(self, 'avec'):
                    ff.write(b'avec 3\n')
                    numpy.savetxt(ff, self.avec)
                if hasattr(self, 'bvec'):
                    ff.write(b'bvec 3\n')
                    numpy.savetxt(ff, self.bvec)

                if hasattr(self, 'avol'):
                    ff.write(b'avol = {0:f}\n'.format(self.avol))
                if hasattr(self, 'bvol'):
                    ff.write(b'bvol = {0:f}\n'.format(self.bvol))
                if hasattr(self, 'natom'):
                    ff.write(b'natom = {0:d}\n'.format(self.natom))
                if hasattr(self, 'nelec'):
                    ff.write(b'nelec = {0:f}\n'.format(self.nelec))

                if hasattr(self, 'rot'):
                    ff.write(b'sym {0:d}\n'.format(self.nsym))
                    numpy.savetxt(ff, self.rot.reshape(self.nsym, 9))

            if level > 1:
                if hasattr(self, 'kunit'):
                    ff.write(b'kunit = {0:s}\n'.format(self.kunit))

                if hasattr(self, 'kpoint'):
                    ff.write(b'kpoint {0:d}\n'.format(self.nkpoint))
                    numpy.savetxt(ff, self.kpoint)

                if hasattr(self, 'kline'):
                    ff.write(b'kline {0:d}\n'.format(self.nkpoint))
                    self.kline.tofile(ff, '\n')

                if hasattr(self, 'kweight'):
                    ff.write(b'kweight {0:d}\n'.format(self.nkpoint))
                    self.kweight.tofile(ff, '\n')

                if hasattr(self, 'kpath'):
                    ff.write(b'kpath {0:d}\n'.format(self.nkpath))
                    numpy.savetxt(ff, self.kpath)

                if hasattr(self, 'kindex'):
                    ff.write(b'kindex {0:d}\n'.format(self.nkpath))
                    self.kindex.tofile(ff, '\n')

                if hasattr(self, 'klabel'):
                    ff.write(b'klabel {0:d}\n'.format(self.nkpath))
                    for ii in range(self.nkpath):
                        ff.write(b'{0:s}\n'.format(self.klabel[ii]))

            if level > 2:
                if hasattr(self, 'eunit'):
                    ff.write(b'eunit = {0:s}\n'.format(self.eunit))

                if hasattr(self, 'energy'):
                    ff.write(s3.format(self.nkpoint, self.nband, self.nspin))
                    self.energy.tofile(ff, '\n')

                if hasattr(self, 'efermi'):
                    ff.write(b'efermi = {0:f}\n'.format(self.efermi))
                if hasattr(self, 'vref'):
                    ff.write(b'vref = {0:f}\n'.format(self.vref))

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

                    ff.write(b'CELL_PARAMETERS {0:s}\n'.format(d1[self.aunit]))
                    numpy.savetxt(ff, self.alat * self.avec)
                    self.set_aunit(oldaunit)

            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpoint'):
                    ff.write(b'K_POINTS {0:s}\n'.format(d2[self.kunit]))
                    ff.write(b'{0:d}\n'.format(self.nkpoint))

                    if hasattr(self, 'kweight'):
                        weight = 2.0 * self.kweight
                    else:
                        weight = numpy.full(self.nkpoint, 2 / self.nkpoint)
                    numpy.savetxt(ff, numpy.concatenate((
                        self.kpoint, weight.reshape(self.nkpoint, 1)), 1))

                if hasattr(self, 'kunit') and hasattr(
                        self, 'kpath') and hasattr(self, 'kindex'):
                    if self.check_kindex() == 'uninitialized':
                        raise ValueError('call calc_kindex')

                    ff.write(b'K_POINTS {0:s}\n'.format(d3[self.kunit]))
                    ff.write(b'{0:d}\n'.format(self.nkpath))

                    knum = numpy.append(numpy.diff(self.kindex), 0)
                    numpy.savetxt(ff, numpy.concatenate((
                        self.kpath, knum.reshape(self.nkpath, 1)), 1))

    def _write_file_wannier_in(self, level, filenames):
        d0 = {'bohr': 'bohr', 'angstrom': 'ang'}
        l0 = [' ', '\n']

        with open(filenames[0], 'wb') as ff:
            if level > 0:
                if hasattr(self, 'aunit') and hasattr(self, 'alat') and hasattr(
                        self, 'avec'):
                    ff.write(b'begin unit_cell_cart\n')
                    ff.write(b'{0:s}'.format(d0[self.aunit]))
                    numpy.savetxt(ff, self.alat * self.avec)
                    ff.write(b'end unit_cell_cart\n')

            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                    oldkunit = self.kunit
                    self.set_kunit('crystal')

                    ff.write(b'begin kpoint_path\n')
                    for ii in range(self.nkpath - 1):
                        for jj in range(2):
                            ff.write(b'{0:s} '.format(self.klabel[ii + jj]))
                            self.kpath[ii + jj, :].tofile(ff, ' ')
                            ff.write(b'{0:s}'.format(l0[jj]))
                    ff.write(b'end kpoint_path\n')

                    self.set_kunit(oldkunit)

    def _write_file_vasp_kpt(self, level, filenames):
        d0 = {'cartesian': 'cartesian', 'crystal': 'reciprocal'}

        with open(filenames[0], 'w') as ff:
            if level > 1:
                if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                    if self.check_kindex() != 'number':
                        raise ValueError("call calc_kindex('number')")

                    ff.write('k-points along high symmetry lines\n')
                    ff.write('{0:d}\n'.format(self.kindex[1]))
                    ff.write('line-mode\n')
                    ff.write('{0:s}'.format(d0[self.kunit]))

                    for ii in range(self.nkpath - 1):
                        ff.write('\n')
                        for jj in range(2):
                            self.kpath[ii + jj, :].tofile(ff, ' ')
                            ff.write(' ! {0:s}\n'.format(self.klabel[ii + jj]))

    def _write_file_lapw_kpt(self, level, filenames, lapwkunit):
        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                if self.check_kindex() == 'uninitialized':
                    raise ValueError('call calc_kindex')

                oldkunit = self.kunit
                self.set_kunit(lapwkunit)
                nline = self.nkpath - 1
                nstep = numpy.diff(self.kindex)
                line = numpy.diff(self.kpath, axis = 0)
                step = line / numpy.broadcast_to(nstep, (3, nline)).transpose()
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

                with open(filenames[0], 'w') as ff:
                    for ii in range(len(labels)):
                        ff.write(s0.format(labels[ii], kpoints[ii], ngrids[
                            ii], weights[ii]))
                    ff.write('END\n')

    def _write_file_boltztrap_in(
            self, level, filenames, deltae, ecut, lpfac, efcut, tmax, deltat,
            ecut2, dosmethod, nband_exclude):
        if level > 2:
            self.set_aunit('bohr')
            self.set_alat(1.0)
            self.set_kunit('crystal')
            self.set_eunit('rydberg')

            with open(filenames[0], 'w') as ff:
                filename_intrans = os.path.basename(filenames[1])
                filename_struct = os.path.basename(filenames[2])
                filename_energy = os.path.basename(filenames[3])

                s0 = ", 'old', 'formatted', 0\n"
                s1 = ", 'unknown', 'formatted', 0\n"
                s2 = ", 'unknown', 'unformatted', 0\n"

                ff.write("5, '{0:s}'{1:s}".format(filename_intrans, s0))
                ff.write("6, '{0:s}.outputtrans'{1:s}".format(self.prefix, s1))
                ff.write("20, '{0:s}'{1:s}".format(filename_struct, s0))
                ff.write("10, '{0:s}'{1:s}".format(filename_energy, s0))
                ff.write("48, '{0:s}.engre'{1:s}".format(self.prefix, s2))
                ff.write("49, '{0:s}.transdos'{1:s}".format(self.prefix, s1))
                ff.write("50, '{0:s}.sigxx'{1:s}".format(self.prefix, s1))
                ff.write("51, '{0:s}.sigxxx'{1:s}".format(self.prefix, s1))
                ff.write("21, '{0:s}.trace'{1:s}".format(self.prefix, s1))
                ff.write("22, '{0:s}.condtens'{1:s}".format(self.prefix, s1))
                ff.write("24, '{0:s}.halltens'{1:s}".format(self.prefix, s1))
                ff.write("30, '{0:s}_BZ.dx'{1:s}".format(self.prefix, s1))
                ff.write("31, '{0:s}_fermi.dx'{1:s}".format(self.prefix, s1))
                ff.write("32, '{0:s}_sigxx.dx'{1:s}".format(self.prefix, s1))
                ff.write("33, '{0:s}_sigyy.dx'{1:s}".format(self.prefix, s1))
                ff.write("34, '{0:s}_sigzz.dx'{1:s}".format(self.prefix, s1))
                ff.write("35, '{0:s}_band.dat'{1:s}".format(self.prefix, s1))
                ff.write("36, '{0:s}_band.gpl'{1:s}".format(self.prefix, s1))
                ff.write("37, '{0:s}_deriv.dat'{1:s}".format(self.prefix, s1))
                ff.write("38, '{0:s}MASS.dat'{1:s}".format(self.prefix, s1))

            with open(filenames[1], 'w') as ff:
                nelec = self.nelec - 2 * nband_exclude

                s0 = '# Format of DOS\n'
                s1 = '# iskip (not presently used) idebug setgap shiftgap\n'
                s2 = '# Fermilevel (Ry), energygrid, energy span around'
                s2 += ' Fermilevel, number of electrons\n'
                s3 = '# CALC (calculate expansion coeff), NOCALC read from'
                s3 += ' file\n'
                s4 = '# lpfac, number of latt-points per k-point\n'
                s5 = '# run mode (only BOLTZ is supported)\n'
                s6 = '# (efcut) energy range of chemical potential\n'
                s7 = '# Tmax, temperature grid\n'
                s8 = '# energyrange of bands given individual DOS output'
                s8 += ' sig_xxx and dos_xxx (xxx is band number)\n'

                ff.write('GENE      {s:0}'.format(s0))
                ff.write('0 0 0 0.0 {s:0}'.format(s1))
                ff.write('{0:f} {1:f} {2:f} {3:f} {4:s}'.format(
                    self.efermi, deltae, ecut, nelec, s2))
                ff.write('CALC {0:s}'.format(s3))
                ff.write('{0:d} {1:s}'.format(lpfac, s4))
                ff.write('BOLTZ {0:s}'.format(s5))
                ff.write('{0:f} {1:s}'.format(efcut, s6))
                ff.write('{0:f} {1:f} {2:s}'.format(tmax, deltat, s7))
                ff.write('{0:f} {1:s}'.format(ecut2, s8))
                ff.write('{0:s}\n'.format(dosmethod))

            with open(filenames[2], 'wb') as ff:
                ff.write(b'{0:s}\n'.format(self.prefix))
                numpy.savetxt(ff, self.avec)
                ff.write(b'{0:d}\n'.format(self.nsym))
                for ir in range(self.nsym):
                    self.rot[ir, :, :].transpose().tofile(ff, ' ')
                    ff.write(b'\n')

            with open(filenames[3], 'w') as ff:
                nband = self.nspin * (self.nband - nband_exclude)
                ff.write('{0:s}\n'.format(self.prefix))
                ff.write('{0:d}\n'.format(self.nkpoint))
                for ii in range(self.nkpoint):
                    self.kpoint[ii, :].tofile(ff, ' ')
                    ff.write(' {0:d}\n'.format(nband))
                    self.energy[ii, nband_exclude:self.nband, :].tofile(
                        ff, '\n')
                    ff.write('\n')

