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
        contents = common._read_file(filenames)

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
                if 'prefix' in line.lower():
                    iprefix = ii
                elif 'aunit' in line.lower():
                    iaunit = ii
                elif 'alat' in line.lower():
                    ialat = ii
                elif 'avec' in line.lower():
                    iavec = ii
                elif 'bvec' in line.lower():
                    ibvec = ii
                elif 'avol' in line.lower():
                    iavol = ii
                elif 'bvol' in line.lower():
                    ibvol = ii
                elif 'natom' in line.lower():
                    inatom = ii
                elif 'nelec' in line.lower():
                    inelec = ii
                elif 'sym' in line.lower():
                    isym = ii
            if level > 1:
                if 'kunit' in line.lower():
                    ikunit = ii
                elif 'kpoint' in line.lower():
                    ikpt = ii
                elif 'kline' in line.lower():
                    ikline = ii
                elif 'kweight' in line.lower():
                    ikweight = ii
                elif 'kpath' in line.lower():
                    ikpath = ii
                elif 'kindex' in line.lower():
                    ikindex = ii
                elif 'klabel' in line.lower():
                    iklabel = ii
            if level > 2:
                if 'eunit' in line.lower():
                    ieunit = ii
                elif 'energy' in line.lower():
                    ienergy = ii
                elif 'efermi' in line.lower():
                    iefrm = ii
                elif 'vref' in line.lower():
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
        contents = common._read_file(filenames)

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
        contents = common._read_file(filenames)

        if level > 1:
            nkpersect = 100
        for ii, line in enumerate(contents[0]):
            if level > 0:
                if 'begin' in line.lower() and 'unit_cell_cart' in line.lower(
                        ):
                    iaunit = ii + 1
            if level > 1:
                if 'bands_num_points' in line.lower():
                    nkpersect = int(line.split()[2])
                elif 'begin' in line.lower() and 'kpoint_path' in line.lower():
                    istart = ii
                elif 'end' in line.lower() and 'kpoint_path' in line.lower():
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
        contents = common._read_file(filenames)

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
        contents = common._read_file(filenames)

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
        contents = common._read_file(filenames)

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
        contents = common._read_file(filenames)

        if level > 1:
            idxkpt = []
        if level > 2:
            idxbnd = []
        ii = 0
        for line in contents[0]:
            if level > 1:
                if ' q =' in line.lower():
                    idxkpt.append(ii)
            if level > 2:
                if 'freq' in line.lower() or 'omega' in line.lower():
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
        contents = common._read_file(filenames)

        if level > 1:
            oldaunit = self.aunit
            self.set_aunit('bohr')
            kscale = self.alat / (2.0 * math.pi)

        if level > 2:
            if etype is None:
                eindex = 6
            else:
                if etype.lower() == 'emf':
                    eindex = 5
                elif etype.lower() == 'eqp':
                    eindex = 6
                else:
                    raise ValueError(etype)

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
        contents = common._read_file(filenames)

        if level > 2:
            if etype is None:
                eindex = 9
            else:
                if etype.lower() == 'edft':
                    eindex = 1
                elif etype.lower() == 'ecor':
                    eindex = 2
                elif etype.lower() == 'eqp0':
                    eindex = 8
                elif etype.lower() == 'eqp1':
                    eindex = 9
                elif etype.lower() == 'eqp0p':
                    eindex = 12
                elif etype.lower() == 'eqp1p':
                    eindex = 13
                else:
                    raise ValueError(etype)

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
        contents = common._read_file(filenames)

        if level > 1:
            oldaunit = self.aunit
            self.set_aunit('angstrom')
            kscale = self.alat / (2.0 * math.pi)

            ndata = len(contents[0])
            for ii in range(ndata):
                if len(contents[0][ii].split()) == 0:
                    nkpoint = ii
                    break

            kline = numpy.empty(nkpoint, float)
            for ikpoint in range(nkpoint):
                words = contents[0][ikpoint].split()
                kline[ikpoint] = float(words[0])

            kline *= kscale
            self.kline = kline
            self.set_aunit(oldaunit)

        if level > 2:
            nband = ndata // (nkpoint + 1)
            nspin = len(contents)

            energy = numpy.empty((nkpoint, nband, nspin), float)
            for ispin in range(nspin):
                ii = 0
                for iband in range(nband):
                    for ikpoint in range(nkpoint):
                        energy[ikpoint, iband, ispin] = float(
                                contents[ispin][ii].split()[1])
                        ii += 1
                    ii += 1

            self.eunit, self.energy = 'ev', energy

    def _write_file_internal(self, level, filenames):
        content = ''

        if level > 0:
            if hasattr(self, 'prefix'):
                content += 'prefix = {0:s}\n'.format(self.prefix)
            if hasattr(self, 'aunit'):
                content += 'aunit = {0:s}\n'.format(self.aunit)
            if hasattr(self, 'alat'):
                content += 'alat = {0:f}\n'.format(self.alat)

            if hasattr(self, 'avec'):
                content += 'avec 3\n'
                for ii in range(3):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(self.avec[
                            ii])

            if hasattr(self, 'bvec'):
                content += 'bvec 3\n'
                for ii in range(3):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(self.bvec[
                            ii])

            if hasattr(self, 'avol'):
                content += 'avol = {0:f}\n'.format(self.avol)
            if hasattr(self, 'bvol'):
                content += 'bvol = {0:f}\n'.format(self.bvol)
            if hasattr(self, 'natom'):
                content += 'natom = {0:d}\n'.format(self.natom)
            if hasattr(self, 'nelec'):
                content += 'nelec = {0:f}\n'.format(self.nelec)

            if hasattr(self, 'rot'):
                content += 'sym {0:d}\n'.format(self.nsym)
                for ii in range(self.nsym):
                    for jj in range(3):
                        content += '{0[0]:2d} {0[1]:2d} {0[2]:2d}'.format(
                                self.rot[ii, jj])
                        if jj < 2:
                            content += '  '
                        else:
                            content += '\n'

        if level > 1:
            if hasattr(self, 'kunit'):
                content += 'kunit = {0:s}\n'.format(self.kunit)

            if hasattr(self, 'kpoint'):
                content += 'kpoint {0:d}\n'.format(self.nkpoint)
                for ikpoint in range(self.nkpoint):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(
                            self.kpoint[ikpoint])

            if hasattr(self, 'kline'):
                content += 'kline {0:d}\n'.format(self.nkpoint)
                for ikpoint in range(self.nkpoint):
                    content += '{0:f}\n'.format(self.kline[ikpoint])

            if hasattr(self, 'kweight'):
                content += 'kweight {0:d}\n'.format(self.nkpoint)
                for ikpoint in range(self.nkpoint):
                    content += '{0:f}\n'.format(self.kweight[ikpoint])

            if hasattr(self, 'kpath'):
                content += 'kpath {0:d}\n'.format(self.nkpath)
                for ikpoint in range(self.nkpath):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(
                            self.kpath[ikpoint])

            if hasattr(self, 'kindex'):
                content += 'kindex {0:d}\n'.format(self.nkpath)
                for ikpoint in range(self.nkpath):
                    content += '{0:n}\n'.format(self.kindex[ikpoint])

            if hasattr(self, 'klabel'):
                content += 'klabel {0:d}\n'.format(self.nkpath)
                for ikpoint in range(self.nkpath):
                    content += '{0:s}\n'.format(self.klabel[ikpoint])

        if level > 2:
            if hasattr(self, 'eunit'):
                content += 'eunit = {0:s}\n'.format(self.eunit)

            if hasattr(self, 'energy'):
                content += 'energy {0:d} {1:d} {2:d}\n'.format(
                        self.nkpoint, self.nband, self.nspin)
                for ikpoint in range(self.nkpoint):
                    for iband in range(self.nband):
                        for ispin in range(self.nspin):
                            content += '{0:f}\n'.format(
                                    self.energy[ikpoint, iband, ispin])

            if hasattr(self, 'efermi'):
                content += 'efermi = {0:f}\n'.format(self.efermi)

            if hasattr(self, 'vref'):
                content += 'vref = {0:f}\n'.format(self.vref)

        common._write_file(filenames, [content])

    def _write_file_pw_in(self, level, filenames):
        content = ''

        if level > 0:
            if hasattr(self, 'aunit') and hasattr(self, 'alat') and hasattr(
                    self, 'avec'):
                content += 'CELL_PARAMETERS '
                if self.aunit == 'bohr':
                    content += 'bohr\n'
                else:
                    content += 'angstrom\n'
                for ii in range(3):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(
                            self.avec[ii] * self.alat)

        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpoint'):
                content += 'K_POINTS '
                if self.kunit == 'cartesian':
                    content += 'tpiba\n'
                else:
                    content += 'crystal\n'
                content += '{0:d}\n'.format(self.nkpoint)
                for ikpoint in range(self.nkpoint):
                    if hasattr(self, 'kweight'):
                        weight = 2.0 * self.kweight[ikpoint]
                    else:
                        weight = 2.0 / float(self.nkpoint)
                    content += '{0[0]:f} {0[1]:f} {0[2]:f} {1:f}\n'.format(
                            self.kpoint[ikpoint], weight)

            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                if self.check_kindex() == 'uninitialized':
                    raise ValueError("call calc_kindex('number'|'density')")

                content += 'K_POINTS '
                if self.kunit == 'cartesian':
                    content += 'tpiba_b\n'
                else:
                    content += 'crystal_b\n'
                content += '{0:d}\n'.format(self.nkpath)
                for ikpoint in range(self.nkpath):
                    if ikpoint < self.nkpath - 1:
                        knum = self.kindex[ikpoint + 1] - self.kindex[ikpoint]
                    else:
                        knum = 0
                    content += '{0[0]:f} {0[1]:f} {0[2]:f} {1:f}\n'.format(
                            self.kpath[ikpoint], float(knum))

        common._write_file(filenames, [content])

    def _write_file_wannier_in(self, level, filenames):
        content = ''

        if level > 0:
            if hasattr(self, 'aunit') and hasattr(self, 'alat') and hasattr(
                    self, 'avec'):
                content += 'begin unit_cell_cart\n'
                if self.aunit == 'bohr':
                    content += 'bohr\n'
                else:
                    content += 'ang\n'
                for ii in range(3):
                    content += '{0[0]:f} {0[1]:f} {0[2]:f}\n'.format(
                            self.avec[ii] * self.alat)
                content += 'end unit_cell_cart\n'

        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                oldkunit = self.kunit
                self.set_kunit('crystal')
                tt = [' ', '\n']

                content += 'begin kpoint_path\n'
                for ikpoint in range(self.nkpath - 1):
                    for ii in range(2):
                        content += ('{0:s} {1[0]:f} {1[1]:f} {1[2]:f}{2:s}'
                                ).format(self.klabel[ikpoint + ii], self.kpath[
                                ikpoint + ii], tt[ii])
                content += 'end kpoint_path\n'

                self.set_kunit(oldkunit)

        common._write_file(filenames, [content])

    def _write_file_vasp_kpt(self, level, filenames):
        content = ''

        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                if self.check_kindex() != 'number':
                    raise ValueError("call calc_kindex('number')")

                content += 'k-points along high symmetry lines\n'
                content += '{0:d}\n'.format(self.kindex[1])
                content += 'line-mode\n'
                if self.kunit == 'cartesian':
                    content += 'cartesian\n'
                else:
                    content += 'reciprocal\n'

                for ikpoint in range(self.nkpath - 1):
                    if ikpoint > 0:
                        content += '\n'
                    for ii in range(2):
                        content += ('{0[0]:f} {0[1]:f} {0[2]:f} ! {1:s}\n'
                                ).format(self.kpath[ikpoint + ii], self.klabel[
                                ikpoint + ii])

        common._write_file(filenames, [content])

    def _write_file_lapw_kpt(self, level, filenames, lapwkunit):
        content = ''

        if level > 1:
            if hasattr(self, 'kunit') and hasattr(self, 'kpath'):
                if self.check_kindex() == 'uninitialized':
                    raise ValueError("call calc_kindex('number'|'density')")

                oldkunit = self.kunit
                self.set_kunit(lapwkunit)

                nline = self.nkpath - 1
                nstep = []
                line = numpy.empty((nline, 3), float)
                step = numpy.empty((nline, 3), float)
                for ii in range(nline):
                    nstep.append(self.kindex[ii + 1] - self.kindex[ii])
                    line[ii] = self.kpath[ii + 1] - self.kpath[ii]
                    step[ii] = line[ii] / float(nstep[ii])
                origin = self.kpath[0] - step[0]

                ngridmax = 8192
                ngrid = []
                for ii in range(nline):
                    denom = []
                    for jj in range(3):
                        frac = fractions.Fraction(step[
                                ii, jj]).limit_denominator()
                        if frac.numerator != 0 and (
                                not frac.denominator in denom):
                            denom.append(frac.denominator)

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
                            label = self.klabel[0]
                        elif jj == nstep[ii]:
                            label = self.klabel[ii + 1]
                        else:
                            label = ''
                        for ll in range(10 - len(label)):
                            label += ' '

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

                for ii in range(len(labels)):
                    content += ('{0:s}{1[0]:5d}{1[1]:5d}{1[2]:5d}{2:5d}{3:s}\n'
                            ).format(labels[ii], kpoints[ii], ngrids[
                            ii], weights[ii])
                content += 'END\n'

                self.set_kunit(oldkunit)

        common._write_file(filenames, [content])

    def _write_file_boltztrap_in(
            self, level, filenames, deltae, ecut, lpfac, efcut, tmax, deltat,
            ecut2, dosmethod, nband_exclude):
        contents = ['', '', '', '']

        if level > 2:
            self.set_aunit('bohr')
            self.set_alat(1.0)
            self.set_kunit('crystal')
            self.set_eunit('rydberg')

            _tail1 = ', \'old\', \'formatted\', 0\n'
            _tail2 = ', \'unknown\', \'formatted\', 0\n'
            _tail3 = ', \'unknown\', \'unformatted\', 0\n'

            filename_intrans = os.path.basename(filenames[1])
            filename_struct = os.path.basename(filenames[2])
            filename_energy = os.path.basename(filenames[3])

            contents[0] += '5, \'' + filename_intrans + '\'' + _tail1
            contents[0] += '6, \'' + self.prefix + '.outputtrans\'' + _tail2
            contents[0] += '20, \'' + filename_struct + '\'' + _tail1
            contents[0] += '10, \'' + filename_energy + '\'' + _tail1
            contents[0] += '48, \'' + self.prefix + '.engre\'' + _tail3
            contents[0] += '49, \'' + self.prefix + '.transdos\'' + _tail2
            contents[0] += '50, \'' + self.prefix + '.sigxx\'' + _tail2
            contents[0] += '51, \'' + self.prefix + '.sigxxx\'' + _tail2
            contents[0] += '21, \'' + self.prefix + '.trace\'' + _tail2
            contents[0] += '22, \'' + self.prefix + '.condtens\'' + _tail2
            contents[0] += '24, \'' + self.prefix + '.halltens\'' + _tail2
            contents[0] += '30, \'' + self.prefix + '_BZ.dx\'' + _tail2
            contents[0] += '31, \'' + self.prefix + '_fermi.dx\'' + _tail2
            contents[0] += '32, \'' + self.prefix + '_sigxx.dx\'' + _tail2
            contents[0] += '33, \'' + self.prefix + '_sigyy.dx\'' + _tail2
            contents[0] += '34, \'' + self.prefix + '_sigzz.dx\'' + _tail2
            contents[0] += '35, \'' + self.prefix + '_band.dat\'' + _tail2
            contents[0] += '36, \'' + self.prefix + '_band.gpl\'' + _tail2
            contents[0] += '37, \'' + self.prefix + '_deriv.dat\'' + _tail2
            contents[0] += '38, \'' + self.prefix + 'MASS.dat\'' + _tail2

            _cmnt1 = ' # Format of DOS\n'
            _cmnt2 = ' # iskip (not presently used) idebug setgap shiftgap\n'
            _cmnt3 = ' # Fermilevel (Ry), energygrid, energy span around'
            _cmnt3 += ' Fermilevel, number of electrons\n'
            _cmnt4 = ' # CALC (calculate expansion coeff), NOCALC read from'
            _cmnt4 += ' file\n'
            _cmnt5 = ' # lpfac, number of latt-points per k-point\n'
            _cmnt6 = ' # run mode (only BOLTZ is supported)\n'
            _cmnt7 = ' # (efcut) energy range of chemical potential\n'
            _cmnt8 = ' # Tmax, temperature grid\n'
            _cmnt9 = ' # energyrange of bands given individual DOS output'
            _cmnt9 += ' sig_xxx and dos_xxx (xxx is band number)\n'
            _cmnt10 = '\n'

            contents[1] += 'GENE     ' + _cmnt1
            contents[1] += '0 0 0 0.0' + _cmnt2
            contents[1] += str(self.efermi) + ' ' + str(deltae) + ' ' + str(
                    ecut) + ' ' + str(self.nelec - float(
                    2 * nband_exclude)) + _cmnt3
            contents[1] += 'CALC' + _cmnt4
            contents[1] += str(lpfac) + _cmnt5
            contents[1] += 'BOLTZ' + _cmnt6
            contents[1] += str(efcut) + _cmnt7
            contents[1] += str(tmax) + ' ' + str(deltat) + _cmnt8
            contents[1] += str(ecut2) + _cmnt9
            contents[1] += dosmethod + _cmnt10

            contents[2] += self.prefix + '\n'
            for ii in range(3):
                for jj in range(3):
                    contents[2] += str(self.avec[ii, jj])
                    if jj < 2:
                        contents[2] += ' '
                    else:
                        contents[2] += '\n'
            contents[2] += str(self.nsym) + '\n'
            for ir in range(self.nsym):
                for ii in range(3):
                    for jj in range(3):
                        contents[2] += str(self.rot[ir, jj, ii])
                        if ii < 2 or jj < 2:
                            contents[2] += ' '
                        else:
                            contents[2] += '\n'

            contents[3] += self.prefix + '\n'
            contents[3] += str(self.nkpoint) + '\n'
            for ikpoint in range(self.nkpoint):
                for ii in range(3):
                    contents[3] += str(self.kpoint[ikpoint, ii]) + ' '
                contents[3] += str(
                        self.nspin * (self.nband - nband_exclude)) + '\n'
                for iband in range(nband_exclude, self.nband):
                    for ispin in range(self.nspin):
                        contents[3] += str(
                                self.energy[ikpoint, iband, ispin]) + '\n'

        common._write_file(filenames, contents)

