"""This module defines class Diagram."""

import numpy as np

import brave.common as common
from brave.energy import Energy
from brave.dos import DOS
from brave.plot import Plot

class Diagram(DOS, Energy):
    """Class for plotting the band diagrams.

    Class Diagram helps to sketch with the proper labeling the electron or
    phonon energy band diagram and/or the electron or phonon density of states.
    """

    @property
    def plot(self):
        """An instance of class Plot holding a sketch of the energy band
    diagram and/or the density of states.
        """
        return self._plot

    @plot.setter
    def plot(self, value):
        if not isinstance(value, Plot):
            raise TypeError('plot {0!r}'.format(value))
        self._plot = value

    @plot.deleter
    def plot(self):
        del self._plot

    def set_plot(
            self, plotformat, fillgap = None, colordos = None, labeldos = None):
        """Sets the new value of plot.

    Args:
        plotformat (str): Plot format. Possible values are 'energy', 'dos' and
            'energy_dos'
        fillgap (list): The band and spin indices of the top valence and bottom
            conduction bands for color filling the band gap. The order is
            [ibandv, ispinv, ibandc, ispinc]. The default is None.
        colordos (list): The colors of PDOS components.
        labeldos (list): The labels of PDOS components.

    You may want to set plot.ylim, plot.note and plot.title before calling
    plot.write.
        """
        _plotformat_list = ['energy', 'dos', 'energy_dos']
        _kind_gap = 'fill'
        _style_gap = ['None', 'None']
        _color_gap = ['yellow', 'none', 'none']
        if not hasattr(self, 'energy') or self.nspin == 1:
            _kind_band = ['plot']
            _style_band = [['solid', 'None']]
            _color_band = [['blue', 'none', 'none']]
            _kind_efermi = 'plot'
            _style_efermi = ['solid', 'None']
            _color_efermi = ['red', 'none', 'none']
            _kind_vref = 'plot'
            _style_vref = ['solid', 'None']
            _color_vref = ['green', 'none', 'none']
        else:
            _kind_band = ['plot', 'plot']
            _style_band = [['solid', 'None'], ['solid', 'None']]
            _color_band = [['red', 'none', 'none'], ['green', 'none', 'none']]
            _kind_efermi = 'plot'
            _style_efermi = ['solid', 'None']
            _color_efermi = ['blue', 'none', 'none']
            _kind_vref = 'plot'
            _style_vref = ['solid', 'None']
            _color_vref = ['black', 'none', 'none']
        _kind_dos = 'plot'
        _style_dos = ['solid', 'None']
        _color_dos = ['blue', 'none', 'none']
        _ypad = 0.1
        _ivlabel_dict = {
                'uc': 'uc$^{-1}$',
                'bohr3': '$a_0^{-3}$',
                'angstrom3': '$\AA^{-3}$',
                'nm3': 'nm$^{-3}$'}
        _ielabel_dict = {
                'ev': 'eV$^{-1}$',
                'rydberg': 'Ry$^{-1}$',
                'hartree': 'Ha$^{-1}$',
                'thz': 'THz$^{-1}$',
                'cm-1': '(cm$^{-1}$)$^{-1}$'}
        _elabel_dict = {
                'ev': 'eV',
                'rydberg': 'Ry',
                'hartree': 'Ha',
                'thz': 'THz',
                'cm-1': 'cm$^{-1}$'}
        _pagesize_list = [[6.0, 4.0], [2.5, 4.0], [7.5, 4.0]]
        _griddim_list = [[1, 1], [1, 1], [1, 5]]
        _gridpos_list = [
                [[0, 1, 0, 1]], [[0, 1, 0, 1]], [[0, 1, 0, 4], [0, 1, 4, 5]]]

        if plotformat in _plotformat_list:
            ii = _plotformat_list.index(plotformat)
            i_e = ii % -2
            i_d = ii - 1
        else:
            raise ValueError(plotformat)

        pagesize = _pagesize_list[ii]
        griddim = _griddim_list[ii]
        gridpos = _gridpos_list[ii]

        data = [[] for ii in range(max(i_e, i_d) + 1)]
        kind = [[] for ii in range(max(i_e, i_d) + 1)]
        style = [[] for ii in range(max(i_e, i_d) + 1)]
        color = [[] for ii in range(max(i_e, i_d) + 1)]
        label = [[] for ii in range(max(i_e, i_d) + 1)]
        xlim = []
        xtick = []
        xgrid = []
        xlabel = []
        ylabel = ['' for ii in range(max(i_e, i_d) + 1)]

        if i_e > -1:
            self.calc_kline()

            if fillgap is not None:
                curve = np.empty((2, self.nkpoint * 2 + 1), float)
                iband = fillgap[0]
                ispin = fillgap[1]
                curve[0, 0:self.nkpoint] = self.kline
                curve[1, 0:self.nkpoint] = self.energy[:, iband, ispin]
                iband = fillgap[2]
                ispin = fillgap[3]
                curve[0, self.nkpoint:2*self.nkpoint] = self.kline[::-1]
                curve[1, self.nkpoint:2*self.nkpoint] = self.energy[
                        ::-1, iband, ispin]
                iband = fillgap[0]
                ispin = fillgap[1]
                curve[0, 2*self.nkpoint] = self.kline[0]
                curve[1, 2*self.nkpoint] = self.energy[0, iband, ispin]
                data[i_e].append(curve)
                kind[i_e].append(_kind_gap)
                style[i_e].append(_style_gap)
                color[i_e].append(_color_gap)
                label[i_e].append('')

            for ispin in range(self.nspin - 1, -1, -1):
                for iband in range(self.nband):
                    curve = np.empty((2, self.nkpoint), float)
                    curve[0, :] = self.kline
                    curve[1, :] = self.energy[:, iband, ispin]
                    data[i_e].append(curve)
                    kind[i_e].append(_kind_band[ispin])
                    style[i_e].append(_style_band[ispin])
                    color[i_e].append(_color_band[ispin])
                    label[i_e].append('')

            if hasattr(self, 'efermi'):
                curve = np.empty((2, 2), float)
                curve[0, 0] = self.kline[0]
                curve[0, 1] = self.kline[-1]
                curve[1, :] = self.efermi
                data[i_e].append(curve)
                kind[i_e].append(_kind_efermi)
                style[i_e].append(_style_efermi)
                color[i_e].append(_color_efermi)
                label[i_e].append('')

            if hasattr(self, 'vref'):
                curve = np.empty((2, 2), float)
                curve[0, 0] = self.kline[0]
                curve[0, 1] = self.kline[-1]
                curve[1, :] = self.vref
                data[i_e].append(curve)
                kind[i_e].append(_kind_vref)
                style[i_e].append(_style_vref)
                color[i_e].append(_color_vref)
                label[i_e].append('')

            xlim.append([0.0, self.kline[self.nkpoint - 1]])

            if hasattr(self, 'kindex') and hasattr(self, 'klabel'):
                _xlist = []
                for ikpoint in range(self.nkpath):
                    _xposition = self.kline[self.kindex[ikpoint]]
                    _xlabel = self.klabel[ikpoint]
                    _xlist.append([_xposition, _xlabel])
                xtick.append(_xlist)
                xgrid.append([['solid', 'black', 1.0], [None, None, None]])
            else:
                xtick.append(None)
                xgrid.append([[None, None, None], [None, None, None]])

            xlabel.append('Wavevector')

        if i_d > -1:
            if i_e > -1:
                self.set_dunit([self.dunit[0], self.eunit])

            ndos = self.dos.shape[0] - 1
            for idos in range(ndos):
                data[i_d].append(self.dos[(idos + 1, 0), :])
                kind[i_d].append(_kind_dos)
                style[i_d].append(_style_dos)
                if colordos is None:
                    color[i_d].append(_color_dos)
                else:
                    color[i_d].append(colordos[idos])
                if labeldos is None:
                    label[i_d].append('')
                else:
                    label[i_d].append(labeldos[idos])

            xlim.append([0.0, np.amax(self.dos[1:, :]) * (1.0 + _ypad)])

            xtick.append(None)
            xgrid.append([[None, None, None], [None, None, None]])

            xlabel.append('DOS ({0:s}{1:s})'.format(_ivlabel_dict[self.dunit[
                    0]], _ielabel_dict[self.dunit[1]]))

        if i_e > -1:
            ylabel[0] = 'Energy ({0:s})'.format(_elabel_dict[self.eunit])
        else:
            ylabel[0] = 'Energy ({0:s})'.format(_elabel_dict[self.dunit[1]])

        plot = Plot()

        plot.data = data
        plot.kind = kind
        plot.style = style
        plot.color = color
        plot.label = label
        plot.xlim = xlim
        plot.xtick = xtick
        plot.xgrid = xgrid
        plot.xlabel = xlabel
        plot.ylabel = ylabel
        plot.pagesize = pagesize
        plot.fontsize = 12.0
        plot.linewidth = 1.0
        plot.labelpad = [2.0, 2.0]
        plot.tickpad = [3.0, 3.0, 3.0, 3.0]
        plot.ticksize = [4.0, 2.0, 4.0, 2.0]
        plot.griddim = griddim
        plot.gridpad = [0.2, 0.0, 0.2]
        plot.gridpos = gridpos

        self.plot = plot

    def read(self, fileformat, filenames, **kwargs):
        """Reads properties from files.

    Args:
        fileformat (str): File format. Possible values are below.
        filenames (list): File names. Possible values are below.

    Inherits fileformat and filenames from classes DOS and Energy.
        """

        super().read(fileformat, filenames, **kwargs)

