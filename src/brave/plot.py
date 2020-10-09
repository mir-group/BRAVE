"""This module defines class Plot."""

import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import brave.common as common

class Plot(object):
    """Class for producing 2D figures."""

    @property
    def data(self):
        """A 2D list of m by n ndarrays of floats holding the data points. The
    first dimension is each plot, the second dimension is each curve, n is the
    number of data points.

    kind        |  m  |  columns of m by n ndarray
    ------------|-----|---------------------------
    'plot'      |  2  |  x and y
    'scatter'   |  3  |  x, y and s
    'bar'       |  6  |  left, height, width, bottom, yerrdn and yerrup
    'errorbar'  |  4  |  x, y, yerr and xerr
    'fill'      |  2  |  x and y
        """
        return self._data

    @data.setter
    def data(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('data {0!r}'.format(value))
        if not all(isinstance(
                item, np.ndarray) and item.dtype == np.dtype(
                'float') and len(item.shape) == 2 and item.shape[0] in (
                2, 3, 4, 6) for group in value for item in group):
            raise ValueError('data {0!r}'.format(value))
        self._data = value

    @data.deleter
    def data(self):
        del self._data

    @property
    def kind(self):
        """A 2D list of strings indicating different types of plots. The first
    dimension is each plot, the second dimension is each curve. Possible values
    are 'plot', 'scatter', 'bar', 'errorbar' and 'fill'. The default is 'plot'.
        """
        return self._kind

    @kind.setter
    def kind(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('kind {0!r}'.format(value))
        if not all(isinstance(item, str) for group in value for item in group):
            raise ValueError('kind {0!r}'.format(value))
        self._kind = value

    @kind.deleter
    def kind(self):
        del self._kind

    @property
    def style(self):
        """A 2D list of length-2 lists of strings holding the line and marker
    styles. The first dimension is each plot, the second dimension is each
    curve, the third dimension is line and marker. Possible values of the line
    style are 'solid', 'dash', 'dot', 'dashdot' and (offset, on-off-dash-seq).
    Possible values of the marker style are 'circle', 'square', 'diamond',
    'triangle_up', 'triangle_down', 'triangle_left', 'triangle_right', 'star',
    'cross' and 'plus'. Each style can be set to any string recognized by
    external library. Each style can be set to 'None' if not used. The default
    is ['solid', 'None'].
        """
        return self._style

    @style.setter
    def style(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('style {0!r}'.format(value))
        if not all(isinstance(item, list) and len(
                item) == 2 for group in value for item in group):
            raise ValueError('style {0!r}'.format(value))
        self._style = value

    @style.deleter
    def style(self):
        del self._style

    @property
    def color(self):
        """A 2D list of length-m lists of strings holding the colors of
    different elements. The first dimension is each plot, the second dimension
    is each curve, the third dimension is a particular element. Possible values
    of the element color are 'black', 'red', 'green', 'blue', 'magenta',
    'yellow', 'cyan' and 'white'. Each color can be set to any string, float or
    tuple recognized by external library. Each color can be set to 'none' if
    not used. The default is ['black', 'none', 'none'].

    kind        |  m  |  elements of length-m list
    ------------|-----|-----------------------------------
    'plot'      |  3  |  line, markeredge and markerface
    'scatter'   |  2  |  markeredge and markerface
    'bar'       |  3  |  barface, baredge and errorbar
    'errorbar'  |  4  |  line, markeredge, markerface and errorbar
    'fill'      |  2  |  patchface and patchedge
        """
        return self._color

    @color.setter
    def color(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('color {0!r}'.format(value))
        if not all(isinstance(item, list) and len(
                item) in (2, 3, 4) for group in value for item in group):
            raise ValueError('color {0!r}'.format(value))
        self._color = value

    @color.deleter
    def color(self):
        del self._color

    @property
    def label(self):
        """A 2D list of strings holding the labels of the curves displayed in
    the plot legends. The first dimension is each plot, the second dimension is
    each curve. A label can be set to an empty string to hide it from the
    legend. The default is None.
        """
        return self._label

    @label.setter
    def label(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('label {0!r}'.format(value))
        if not all(isinstance(item, str) or (
                item is None) for group in value for item in group):
            raise ValueError('label {0!r}'.format(value))
        self._label = value

    @label.deleter
    def label(self):
        del self._label

    @property
    def zorder(self):
        """A 2D list of integers holding the drawing orders of the curves. The
    first dimension is each plot, the second dimension is each curve. The
    curves with lower zorder values are drawn first. The default is 2.
        """
        return self._zorder

    @zorder.setter
    def zorder(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('zorder {0!r}'.format(value))
        if not all(isinstance(item, int) for group in value for item in group):
            raise ValueError('zorder {0!r}'.format(value))
        self._zorder = value

    @zorder.deleter
    def zorder(self):
        del self._zorder

    @property
    def legend(self):
        """A list of [string, integer, [float, float, float, float]] specifying
    the parameters of the plot legends. The first dimension is each plot, the
    string is the location of the legend, the integer is the number of columns
    in the legend, the floats are the left and bottom positions and the width
    and height of the legend. Can be set to None to hide the legend. The
    default is ['best', 1, None].
        """
        return self._legend

    @legend.setter
    def legend(self, value):
        if not isinstance(value, list):
            raise TypeError('legend {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 3) or item is None for item in value):
            raise ValueError('legend {0!r}'.format(value))
        self._legend = value

    @legend.deleter
    def legend(self):
        del self._legend

    @property
    def xscale(self):
        """A list of strings holding the scalings of the x-axes. The first
    dimension is each plot. The possible values are 'linear' and 'log'. Each
    scaling can be set to any string recongized by external library. Can be set
    to None to use the default. The default is 'linear'.
        """
        return self._xscale

    @xscale.setter
    def xscale(self, value):
        if not isinstance(value, list):
            raise TypeError('xscale {0!r}'.format(value))
        if not all(isinstance(item, str) or item is None for item in value):
            raise ValueError('xscale {0!r}'.format(value))
        self._xscale = value

    @xscale.deleter
    def xscale(self):
        del self._xscale

    @property
    def yscale(self):
        """A list of strings holding the scalings of the y-axes. The first
    dimension is each plot. The possible values are 'linear' and 'log'. Each
    scaling can be set to any string recongized by external library. Can be set
    to None to use the default. The default is 'linear'.
        """
        return self._yscale

    @yscale.setter
    def yscale(self, value):
        if not isinstance(value, list):
            raise TypeError('yscale {0!r}'.format(value))
        if not all(isinstance(item, str) or item is None for item in value):
            raise ValueError('yscale {0!r}'.format(value))
        self._yscale = value

    @yscale.deleter
    def yscale(self):
        del self._yscale

    @property
    def xlim(self):
        """A list of length-2 lists of floats holding the limits of the x-axes.
    The first dimension is each plot, the floats are the starting and ending
    values on the x-axis. Can be set to None to autoscale the x-axis. The
    default is None.
        """
        return self._xlim

    @xlim.setter
    def xlim(self, value):
        if not isinstance(value, list):
            raise TypeError('xlim {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('xlim {0!r}'.format(value))
        self._xlim = value

    @xlim.deleter
    def xlim(self):
        del self._xlim

    @property
    def ylim(self):
        """A list of length-2 lists of floats holding the limits of the y-axes.
    The first dimension is each plot, the floats are the starting and ending
    values on the y-axis. Can be set to None to autoscale the y-axis. The
    default is None.
        """
        return self._ylim

    @ylim.setter
    def ylim(self, value):
        if not isinstance(value, list):
            raise TypeError('ylim {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('ylim {0!r}'.format(value))
        self._ylim = value

    @ylim.deleter
    def ylim(self):
        del self._ylim

    @property
    def xdel(self):
        """A list of length-2 lists of floats holding the intervals for the
    ticks on the x-axes. The first dimension is each plot, the floats are the
    intervals for the major and minor ticks on the x-axis. Can be set to None
    to use automatic ticks on the x-axis. The default is None.
        """
        return self._xdel

    @xdel.setter
    def xdel(self, value):
        if not isinstance(value, list):
            raise TypeError('xdel {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('xdel {0!r}'.format(value))
        self._xdel = value

    @xdel.deleter
    def xdel(self):
        del self._xdel

    @property
    def ydel(self):
        """A list of length-2 lists of floats holding the intervals for the
    ticks on the y-axes. The first dimension is each plot, the floats are the
    intervals for the major and minor ticks on the y-axis. Can be set to None
    to use automatic ticks on the y-axis. The default is None.
        """
        return self._ydel

    @ydel.setter
    def ydel(self, value):
        if not isinstance(value, list):
            raise TypeError('ydel {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('ydel {0!r}'.format(value))
        self._ydel = value

    @ydel.deleter
    def ydel(self):
        del self._ydel

    @property
    def xtick(self):
        """A 2D list of [float, string] holding the custom ticks on the x-axes.
    The first dimension is each plot, the second dimension is each tick on the
    x-axis, the float and the string are the position and the label of the
    tick. Set to None to use automatic ticks on the x-axis. The default is
    None.
        """
        return self._xtick

    @xtick.setter
    def xtick(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) or item is None for item in value):
            raise TypeError('xtick {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 2 for group in (
                value) if group is not None for item in group):
            raise ValueError('xtick {0!r}'.format(value))
        self._xtick = value

    @xtick.deleter
    def xtick(self):
        del self._xtick

    @property
    def ytick(self):
        """A 2D list of [float, string] holding the custom ticks on the y-axes.
    The first dimension is each plot, the second dimension is each tick on the
    y-axis, the float and the string are the position and the label of the
    tick. Set to None to use automatic ticks on the y-axis. The default is
    None.
        """
        return self._ytick

    @ytick.setter
    def ytick(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) or item is None for item in value):
            raise TypeError('ytick {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 2 for group in (
                value) if group is not None for item in group):
            raise ValueError('ytick {0!r}'.format(value))
        self._ytick = value

    @ytick.deleter
    def ytick(self):
        del self._ytick

    @property
    def xgrid(self):
        """A 2D list of [string, string, float] holding the vertical grid
    lines. The first dimension is each plot, the second dimension is the major
    and minor grids, the strings and the float are the style, the color, and
    the width of the grid lines. Set to None to disable the vertical grid
    lines. The default is None.
        """
        return self._xgrid

    @xgrid.setter
    def xgrid(self, value):
        if not isinstance(value, list):
            raise TypeError('xgrid {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 2 and all((
                isinstance(elem, list) and len(elem) == 3) or elem is (
                None) for elem in item) for item in value):
            raise ValueError('xgrid {0!r}'.format(value))
        self._xgrid = value

    @xgrid.deleter
    def xgrid(self):
        del self._xgrid

    @property
    def ygrid(self):
        """A 2D list of [string, string, float] holding the horizontal grid
    lines. The first dimension is each plot, the second dimension is the major
    and minor grids, the strings and the float are the style, the color, and
    the width of the grid lines. Set to None to disable the horizontal grid
    lines. The default is None.
        """
        return self._ygrid

    @ygrid.setter
    def ygrid(self, value):
        if not isinstance(value, list):
            raise TypeError('ygrid {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 2 and all((
                isinstance(elem, list) and len(elem) == 3) or elem is (
                None) for elem in item) for item in value):
            raise ValueError('ygrid {0!r}'.format(value))
        self._ygrid = value

    @ygrid.deleter
    def ygrid(self):
        del self._ygrid

    @property
    def xlabel(self):
        """A list of strings holding the labels of the x-axes. The first
    dimension is each plot. A label can be set to an empty string to hide it
    from the plot. The default is None.
        """
        return self._xlabel

    @xlabel.setter
    def xlabel(self, value):
        if not isinstance(value, list):
            raise TypeError('xlabel {0!r}'.format(value))
        if not all(isinstance(item, str) for item in value):
            raise ValueError('xlabel {0!r}'.format(value))
        self._xlabel = value

    @xlabel.deleter
    def xlabel(self):
        del self._xlabel

    @property
    def ylabel(self):
        """A list of strings holding the labels of the y-axes. The first
    dimension is each plot. A label can be set to an empty string to hide it
    from the plot. The default is None.
        """
        return self._ylabel

    @ylabel.setter
    def ylabel(self, value):
        if not isinstance(value, list):
            raise TypeError('ylabel {0!r}'.format(value))
        if not all(isinstance(item, str) for item in value):
            raise ValueError('ylabel {0!r}'.format(value))
        self._ylabel = value

    @ylabel.deleter
    def ylabel(self):
        del self._ylabel

    @property
    def xlabpos(self):
        """A list of [float, string] holding the positions of the labels for
    the x-axes. The first dimension is each plot, the float and the string are
    the locations and the alignments of the label for the x-axis. Possible
    values of the location are from 0 to 1, of the alignment are 'center',
    'left' and 'right'. Set to None to automatically place the label on the
    x-axis. The default is None.
        """
        return self._xlabpos

    @xlabpos.setter
    def xlabpos(self, value):
        if not isinstance(value, list):
            raise TypeError('xlabpos {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('xlabpos {0!r}'.format(value))
        self._xlabpos = value

    @xlabpos.deleter
    def xlabpos(self):
        del self._xlabpos

    @property
    def ylabpos(self):
        """A list of [float, string] holding the positions of the labels for
    the y-axes. The first dimension is each plot, the float and the string are
    the locations and the alignments of the label for the y-axis. Possible
    values of the location are from 0 to 1, of the alignment are 'center',
    'left' and 'right'. Set to None to automatically place the label on the
    y-axis. The default is None.
        """
        return self._ylabpos

    @ylabpos.setter
    def ylabpos(self, value):
        if not isinstance(value, list):
            raise TypeError('ylabpos {0!r}'.format(value))
        if not all((isinstance(item, list) and len(
                item) == 2) or item is None for item in value):
            raise ValueError('ylabpos {0!r}'.format(value))
        self._ylabpos = value

    @ylabpos.deleter
    def ylabpos(self):
        del self._ylabpos

    @property
    def note(self):
        """A 2D list of [float, float, string, string, string, string, float]
    holding the custom labels on the plots. The first dimension is each plot,
    the second dimension is each label, the floats are the locations along the
    x-axis and the y-axis and the font size scaling factor, the strings are the
    horizontal and vertical alignments, the text and the color. Locations
    (0, 0) and (1, 1) correspond to the bottom left and the top right corners
    of the plot. Possible values of the horizontal alignment are 'left' and
    'right', of the vertical alignment are 'bottom' and 'top'. Set to an empty
    list to hide it from the plot. The default is None.
        """
        return self._note

    @note.setter
    def note(self, value):
        if not isinstance(value, list) or not all(isinstance(
                item, list) for item in value):
            raise TypeError('note {0!r}'.format(value))
        if not all(isinstance(item, list) and len(
                item) == 7 for group in value for item in group):
            raise ValueError('note {0!r}'.format(value))
        self._note = value

    @note.deleter
    def note(self):
        del self._note

    @property
    def title(self):
        """A list of [float, float, string, string, string, string, float]
    holding the custom labels on the figure. The first dimension is each label,
    the floats are the locations along the x-axis and the y-axis and the font
    size scaling factor, the strings are the horizontal and vertical
    alignments, the text and the color. Locations (0, 0) and (1, 1) correspond
    to the bottom left and the top right corners of the figure. Possible values
    of the horizontal alignment are 'left' and 'right', of the vertical
    alignment are 'bottom' and 'top'. Set to an empty list to hide it from the
    figure. The default is None.
        """
        return self._title

    @title.setter
    def title(self, value):
        if not isinstance(value, list):
            raise TypeError('title {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 7 for item in value):
            raise ValueError('title {0!r}'.format(value))
        self._title = value

    @title.deleter
    def title(self):
        del self._title

    @property
    def pagesize(self):
        """A length-2 list of floats holding the width and the height of the
    figure in inches. The default is [8.0, 6.0].
        """
        return self._pagesize

    @pagesize.setter
    def pagesize(self, value):
        if not isinstance(value, list):
            raise TypeError('pagesize {0!r}'.format(value))
        if len(value) != 2 or not all(isinstance(
                item, float) for item in value):
            raise ValueError('pagesize {0!r}'.format(value))
        self._pagesize = value

    @pagesize.deleter
    def pagesize(self):
        del self._pagesize

    @property
    def fontsize(self):
        """A float holding the size of the font in points. The default is 18.0.
        """
        return self._fontsize

    @fontsize.setter
    def fontsize(self, value):
        if not isinstance(value, float):
            raise TypeError('fontsize {0!r}'.format(value))
        self._fontsize = value

    @fontsize.deleter
    def fontsize(self):
        del self._fontsize

    @property
    def linewidth(self):
        """A float holding the width of the line in points. The default is 1.5.
        """
        return self._linewidth

    @linewidth.setter
    def linewidth(self, value):
        if not isinstance(value, float):
            raise TypeError('linewidth {0!r}'.format(value))
        self._linewidth = value

    @linewidth.deleter
    def linewidth(self):
        del self._linewidth

    @property
    def markersize(self):
        """A float holding the size of the marker in points. The default is 8.0.
        """
        return self._markersize

    @markersize.setter
    def markersize(self, value):
        if not isinstance(value, float):
            raise TypeError('markersize {0!r}'.format(value))
        self._markersize = value

    @markersize.deleter
    def markersize(self):
        del self._markersize

    @property
    def labelpad(self):
        """A length-2 list of floats holding the amount of padding for the
    labels on the x-axis and the y-axis in points. The default is [6.0, 6.0].
        """
        return self._labelpad

    @labelpad.setter
    def labelpad(self, value):
        if not isinstance(value, list):
            raise TypeError('labelpad {0!r}'.format(value))
        if len(value) != 2 or not all(isinstance(
                item, float) for item in value):
            raise ValueError('labelpad {0!r}'.format(value))
        self._labelpad = value

    @labelpad.deleter
    def labelpad(self):
        del self._labelpad

    @property
    def tickpad(self):
        """A length-4 list of floats holding the amount of padding for the
    major and minor ticks on the x-axis and the major and minor ticks on the
    y-axis in points. The default is [6.0, 6.0, 6.0, 6.0].
        """
        return self._tickpad

    @tickpad.setter
    def tickpad(self, value):
        if not isinstance(value, list):
            raise TypeError('tickpad {0!r}'.format(value))
        if len(value) != 4 or not all(isinstance(
                item, float) for item in value):
            raise ValueError('tickpad {0!r}'.format(value))
        self._tickpad = value

    @tickpad.deleter
    def tickpad(self):
        del self._tickpad

    @property
    def ticksize(self):
        """A length-4 list of floats holding the size of the major and minor
    ticks on the x-axis and the major and minor ticks on the y-axis in points.
    The default is [6.0, 3.0, 6.0, 3.0].
        """
        return self._ticksize

    @ticksize.setter
    def ticksize(self, value):
        if not isinstance(value, list):
            raise TypeError('ticksize {0!r}'.format(value))
        if len(value) != 4 or not all(isinstance(
                item, float) for item in value):
            raise ValueError('ticksize {0!r}'.format(value))
        self._ticksize = value

    @ticksize.deleter
    def ticksize(self):
        del self._ticksize

    @property
    def griddim(self):
        """A length-2 list of integers holding the number of rows and columns
    in the grid whose cells hold the plots. The cells are numbered from the top
    left corner starting from zero. The default is [1, 1].
        """
        return self._griddim

    @griddim.setter
    def griddim(self, value):
        if not isinstance(value, list):
            raise TypeError('griddim {0!r}'.format(value))
        if len(value) != 2 or not all(isinstance(item, int) for item in value):
            raise ValueError('griddim {0!r}'.format(value))
        self._griddim = value

    @griddim.deleter
    def griddim(self):
        del self._griddim

    @property
    def gridpad(self):
        """A length-3 list of floats holding the padding of the figure and the
    vertical and horizontal padding of the plots in units of fontsize. The
    default is [0.5, 0.0, 0.0].
    """
        return self._gridpad

    @gridpad.setter
    def gridpad(self, value):
        if not isinstance(value, list):
            raise TypeError('gridpad {0!r}'.format(value))
        if len(value) != 3 or not all(isinstance(
                item, float) for item in value):
            raise ValueError('gridpad {0!r}'.format(value))
        self._gridpad = value

    @gridpad.deleter
    def gridpad(self):
        del self._gridpad

    @property
    def gridpos(self):
        """A list of length-4 lists of integers holding the positions of the
    plots in the grid. The first dimension is each plot, the integers are the
    top, bottom, left and right borders of the plot in the grid. The default is
    [[0, 1, 0, 1]].
        """
        return self._gridpos

    @gridpos.setter
    def gridpos(self, value):
        if not isinstance(value, list):
            raise TypeError('gridpos {0!r}'.format(value))
        if not all(isinstance(item, list) and len(item) == 4 and all(
                isinstance(elem, int) for elem in item) for item in value):
            raise ValueError('gridpos {0!r}'.format(value))
        self._gridpos = value

    @gridpos.deleter
    def gridpos(self):
        del self._gridpos

    @property
    def frame(self):
        """A boolean indicating whether to display the legend frame. The
    default is False.
        """
        return self._frame

    @frame.setter
    def frame(self, value):
        if not isinstance(value, bool):
            raise TypeError('frame {0!r}'.format(value))
        self._frame = value

    @frame.deleter
    def frame(self):
        del self._frame

    def write(self, libraryname, fileformat, filename):
        """Writes the plot to the file.

    Args:
        libraryname (str): Name of the plotting library. Currently the only
            possible value is 'matplotlib'.
        fileformat (list): File format. Possible values are 'png', 'eps' and
            'pdf'.
        filename (str): File name.
        """

        if libraryname == 'matplotlib':
            self._writep_matplotlib(fileformat, filename)
        else:
            raise ValueError(libraryname)

    def _writep_matplotlib(self, fileformat, filename):

        style_dict = [{
            'solid': '-',
            'dash': '--',
            'dot': ':',
            'dashdot': '-.'}, {
            'circle': 'o',
            'square': 's',
            'diamond': 'D',
            'triangle_up': '^',
            'triangle_down': 'v',
            'triangle_left': '<',
            'triangle_right': '>',
            'star': '*',
            'cross': 'x',
            'plus': '+'}]

        color_dict = {
            'black': 'k',
            'red': 'r',
            'green': 'g',
            'blue': 'b',
            'magenta': 'm',
            'yellow': 'y',
            'cyan': 'c',
            'white': 'w'}

        matplotlib.rcParams['lines.linewidth'] = self.linewidth
        matplotlib.rcParams['lines.markeredgewidth'] = self.linewidth
        matplotlib.rcParams['lines.markersize'] = self.markersize
        matplotlib.rcParams['font.size'] = self.fontsize
        matplotlib.rcParams['axes.linewidth'] = self.linewidth
        matplotlib.rcParams['axes.titlesize'] = self.fontsize
        matplotlib.rcParams['axes.formatter.limits'] = -4, 4
        matplotlib.rcParams['xtick.top'] = True
        matplotlib.rcParams['xtick.bottom'] = True
        matplotlib.rcParams['xtick.major.size'] = self.ticksize[0]
        matplotlib.rcParams['xtick.minor.size'] = self.ticksize[1]
        matplotlib.rcParams['xtick.major.width'] = self.linewidth
        matplotlib.rcParams['xtick.minor.width'] = self.linewidth
        matplotlib.rcParams['xtick.major.pad'] = self.tickpad[0]
        matplotlib.rcParams['xtick.minor.pad'] = self.tickpad[1]
        matplotlib.rcParams['xtick.direction'] = 'in'
        matplotlib.rcParams['ytick.left'] = True
        matplotlib.rcParams['ytick.right'] = True
        matplotlib.rcParams['ytick.major.size'] = self.ticksize[2]
        matplotlib.rcParams['ytick.minor.size'] = self.ticksize[3]
        matplotlib.rcParams['ytick.major.width'] = self.linewidth
        matplotlib.rcParams['ytick.minor.width'] = self.linewidth
        matplotlib.rcParams['ytick.major.pad'] = self.tickpad[2]
        matplotlib.rcParams['ytick.minor.pad'] = self.tickpad[3]
        matplotlib.rcParams['ytick.direction'] = 'in'
        matplotlib.rcParams['grid.linewidth'] = self.linewidth
        matplotlib.rcParams['errorbar.capsize'] = 3.0
        matplotlib.rcParams['legend.framealpha'] = 1.0
        matplotlib.rcParams['legend.edgecolor'] = 'k'
        matplotlib.rcParams['legend.fancybox'] = False
        matplotlib.rcParams['legend.numpoints'] = 1
        matplotlib.rcParams['legend.scatterpoints'] = 1
        matplotlib.rcParams['legend.markerscale'] = 1.0
        matplotlib.rcParams['legend.fontsize'] = self.fontsize
        if self.frame:
            matplotlib.rcParams['legend.frameon'] = True
            matplotlib.rcParams['legend.borderpad'] = 0.5
        else:
            matplotlib.rcParams['legend.frameon'] = False
            matplotlib.rcParams['legend.borderpad'] = 0.2
        matplotlib.rcParams['legend.labelspacing'] = 0.3
        matplotlib.rcParams['legend.handlelength'] = 1.0
        matplotlib.rcParams['legend.handleheight'] = 0.5
        matplotlib.rcParams['legend.handletextpad'] = 0.5
        matplotlib.rcParams['legend.columnspacing'] = 1.0
        matplotlib.rcParams['savefig.dpi'] = 300.0
        matplotlib.rcParams['savefig.pad_inches'] = self.gridpad[
                0] * self.fontsize / 72.0

        fig1 = plt.figure(figsize = (self.pagesize[0], self.pagesize[1]))
        gs1 = gridspec.GridSpec(self.griddim[0], self.griddim[1])

        artist1 = []
        for ii, data in enumerate(self.data):
            ax1 = fig1.add_subplot(
                gs1[self.gridpos[ii][0]:self.gridpos[ii][1],
                self.gridpos[ii][2]:self.gridpos[ii][3]], axisbelow = True)

            legend1 = False
            for jj, curve in enumerate(data):
                if hasattr(self, 'kind'):
                    kind1 = self.kind[ii][jj]
                else:
                    kind1 = 'plot'

                if hasattr(self, 'style'):
                    style1 = self.style[ii][jj]
                else:
                    style1 = ['solid', 'None']

                if hasattr(self, 'color'):
                    color1 = self.color[ii][jj]
                else:
                    color1 = ['black', 'none', 'none']

                if hasattr(self, 'label'):
                    label1 = self.label[ii][jj]
                    if label1 != '' and label1 is not None:
                        legend1 = True
                else:
                    label1 = ''

                if hasattr(self, 'zorder'):
                    zorder1 = self.zorder[ii][jj]
                else:
                    zorder1 = 2

                style2 = []
                for kk in range(2):
                    aa = style1[kk]
                    bb = aa
                    if isinstance(aa, str):
                        if aa in style_dict[kk]:
                            bb = style_dict[kk][aa]
                    style2.append(bb)

                color2 = []
                for cc in color1:
                    dd = cc
                    if isinstance(cc, str):
                        if cc in color_dict:
                            dd = color_dict[cc]
                    color2.append(dd)

                if kind1 == 'plot':
                    curve1 = ax1.plot(
                        curve[0, :], curve[1, :], linestyle = style2[0],
                        marker = style2[1], color = color2[0],
                        markeredgecolor = color2[1], markerfacecolor = color2[
                        2], label = label1, zorder = zorder1)

                elif kind1 == 'scatter':
                    curve1 = ax1.scatter(
                        curve[0, :], curve[1, :], s = curve[2, :],
                        marker = style2[1], edgecolors = color2[0],
                        color = color2[1], label = label1, zorder = zorder1)

                elif kind1 == 'bar':
                    if max(np.amax(curve[4, :]), np.amax(
                            curve[5, :])) > common.EPS12:
                        yerr = [curve[4, :], curve[5, :]]
                    else:
                        yerr = None
                    curve1 = ax1.bar(
                        curve[0, :], curve[1, :], width = curve[2, :],
                        bottom = curve[3, :], yerr = yerr, color = color2[0],
                        edgecolor = color2[1], ecolor = color2[2],
                        linewidth = self.linewidth, label = label1,
                        zorder = zorder1)

                elif kind1 == 'errorbar':
                    if np.amax(curve[2, :]) > common.EPS12:
                        yerr = curve[2, :]
                    else:
                        yerr = None
                    if np.amax(curve[3, :]) > common.EPS12:
                        xerr = curve[3, :]
                    else:
                        xerr = None
                    curve1 = ax1.errorbar(
                        curve[0, :], curve[1, :], yerr = yerr, xerr = xerr,
                        linestyle = style2[0], marker = style2[1],
                        color = color2[0], markeredgecolor = color2[1],
                        markerfacecolor = color2[2], ecolor = color2[3],
                        label = label1, zorder = zorder1)

                elif kind1 == 'fill':
                    curve1 = ax1.fill(
                        curve[0, :], curve[1, :], linestyle = style2[0],
                        facecolor = color2[0], edgecolor = color2[1],
                        label = label1, zorder = zorder1)

                else:
                    raise ValueError('kind {0!r}'.format(kind1))

            if legend1:
                if hasattr(self, 'legend'):
                    if self.legend[ii] is not None:
                        if isinstance(self.legend[ii][0], str):
                            loc1 = self.legend[ii][0]
                        else:
                            loc1 = 'best'
                        if isinstance(self.legend[ii][1], int):
                            ncol1 = self.legend[ii][1]
                        else:
                            ncol1 = 1
                        if isinstance(self.legend[ii][2], list):
                            pad1 = 0.0
                            ancor1 = self.legend[ii][2]
                            tran1 = fig1.transFigure
                        else:
                            pad1 = 0.5
                            ancor1 = None
                            tran1 = ax1.transAxes
                    else:
                        loc1 = 'best'
                        ncol1 = 1
                        pad1 = 0.5
                        ancor1 = None
                        tran1 = ax1.transAxes
                else:
                    loc1 = 'best'
                    ncol1 = 1
                    pad1 = 0.5
                    ancor1 = None
                    tran1 = ax1.transAxes
                legend2 = ax1.legend(
                        loc = loc1, ncol = ncol1, borderaxespad = pad1,
                        bbox_to_anchor = ancor1, bbox_transform = tran1)
                if self.frame:
                    legend2.get_frame().set_linewidth(self.linewidth)

            if hasattr(self, 'xscale'):
                ax1.set_xscale(self.xscale[ii])

            if hasattr(self, 'yscale'):
                ax1.set_yscale(self.yscale[ii])

            ax1.autoscale(enable = True, axis = 'both', tight = False)

            if hasattr(self, 'xlim'):
                if self.xlim[ii] is not None:
                    ax1.set_xlim(self.xlim[ii][0], self.xlim[ii][1])

            if hasattr(self, 'ylim'):
                if self.ylim[ii] is not None:
                    ax1.set_ylim(self.ylim[ii][0], self.ylim[ii][1])

            xmin, xmax, ymin, ymax = ax1.axis()

            if hasattr(self, 'xdel'):
                if self.xdel[ii] is not None:
                    ax1.xaxis.set_major_locator(
                        matplotlib.ticker.MultipleLocator(self.xdel[ii][0]))
                    ax1.xaxis.set_minor_locator(
                        matplotlib.ticker.MultipleLocator(self.xdel[ii][1]))

            if hasattr(self, 'ydel'):
                if self.ydel[ii] is not None:
                    ax1.yaxis.set_major_locator(
                        matplotlib.ticker.MultipleLocator(self.ydel[ii][0]))
                    ax1.yaxis.set_minor_locator(
                        matplotlib.ticker.MultipleLocator(self.ydel[ii][1]))

            if hasattr(self, 'xtick'):
                if self.xtick[ii] is not None:
                    position2 = []
                    label2 = []
                    for item in self.xtick[ii]:
                        position2.append(item[0])
                        label2.append(item[1])
                    ax1.set_xticks(position2)
                    ax1.set_xticklabels(label2)

            if hasattr(self, 'ytick'):
                if self.ytick[ii] is not None:
                    position2 = []
                    label2 = []
                    for item in self.ytick[ii]:
                        position2.append(item[0])
                        label2.append(item[1])
                    ax1.set_yticks(position2)
                    ax1.set_yticklabels(label2)

            if hasattr(self, 'xgrid'):
                if self.xgrid[ii][0] is not None:
                    style1 = self.xgrid[ii][0][0]
                    color1 = self.xgrid[ii][0][1]
                    width2 = self.xgrid[ii][0][2]
                    if not None in (style1, color1, width2):
                        style2 = style1
                        if isinstance(style1, str):
                            if style1 in style_dict[0]:
                                style2 = style_dict[0][style1]
                        color2 = color1
                        if isinstance(color1, str):
                            if color1 in color_dict:
                                color2 = color_dict[color1]
                        ax1.grid(
                            b = True, which = 'major', axis = 'x',
                            linestyle = style2, color = color2,
                            linewidth = width2)

                if self.xgrid[ii][1] is not None:
                    style1 = self.xgrid[ii][1][0]
                    color1 = self.xgrid[ii][1][1]
                    width2 = self.xgrid[ii][1][2]
                    if not None in (style1, color1, width2):
                        style2 = style1
                        if isinstance(style1, str):
                            if style1 in style_dict[0]:
                                style2 = style_dict[0][style1]
                        color2 = color1
                        if isinstance(color1, str):
                            if color1 in color_dict:
                                color2 = color_dict[color1]
                        ax1.grid(
                            b = True, which = 'minor', axis = 'x',
                            linestyle = style2, color = color2,
                            linewidth = width2)

            if hasattr(self, 'ygrid'):
                if self.ygrid[ii][0] is not None:
                    style1 = self.ygrid[ii][0][0]
                    color1 = self.ygrid[ii][0][1]
                    width2 = self.ygrid[ii][0][2]
                    if not None in (style1, color1, width2):
                        style2 = style1
                        if isinstance(style1, str):
                            if style1 in style_dict[0]:
                                style2 = style_dict[0][style1]
                        color2 = color1
                        if isinstance(color1, str):
                            if color1 in color_dict:
                                color2 = color_dict[color1]
                        ax1.grid(
                            b = True, which = 'major', axis = 'y',
                            linestyle = style2, color = color2,
                            linewidth = width2)

                if self.ygrid[ii][1] is not None:
                    style1 = self.ygrid[ii][1][0]
                    color1 = self.ygrid[ii][1][1]
                    width2 = self.ygrid[ii][1][2]
                    if not None in (style1, color1, width2):
                        style2 = style1
                        if isinstance(style1, str):
                            if style1 in style_dict[0]:
                                style2 = style_dict[0][style1]
                        color2 = color1
                        if isinstance(color1, str):
                            if color1 in color_dict:
                                color2 = color_dict[color1]
                        ax1.grid(
                            b = True, which = 'minor', axis = 'y',
                            linestyle = style2, color = color2,
                            linewidth = width2)

            if hasattr(self, 'xlabel'):
                if self.xlabel[ii] == '':
                    label2 = []
                    for label3 in ax1.get_xticklabels():
                        label2.append('')
                    ax1.set_xticklabels(label2)
                else:
                    position2 = 0.5
                    align2 = 'center'
                    if hasattr(self, 'xlabpos'):
                        if self.xlabpos[ii] is not None:
                            position2 = self.xlabpos[ii][0]
                            align2 = self.xlabpos[ii][1]
                    label3 = ax1.set_xlabel(
                        self.xlabel[ii], labelpad = self.labelpad[0],
                        multialignment = 'center', x = position2, ha = align2)
                    artist1.append(label3)

            if hasattr(self, 'ylabel'):
                if self.ylabel[ii] == '':
                    label2 = []
                    for label3 in ax1.get_yticklabels():
                        label2.append('')
                    ax1.set_yticklabels(label2)
                else:
                    position2 = 0.5
                    align2 = 'center'
                    if hasattr(self, 'ylabpos'):
                        if self.ylabpos[ii] is not None:
                            position2 = self.ylabpos[ii][0]
                            align2 = self.ylabpos[ii][1]
                    label3 = ax1.set_ylabel(
                        self.ylabel[ii], labelpad = self.labelpad[1],
                        multialignment = 'center', y = position2, ha = align2)
                    artist1.append(label3)

            if hasattr(self, 'note'):
                for item in self.note[ii]:
                    position1 = item[0]
                    position2 = item[1]
                    align1 = item[2]
                    align2 = item[3]
                    text1 = item[4]
                    color1 = item[5]
                    fontscale1 = item[6]
                    text2 = ax1.text(
                        position1, position2, text1,
                        horizontalalignment = align1,
                        verticalalignment = align2,
                        transform = ax1.transAxes, color = color1,
                        fontsize = fontscale1 * self.fontsize)
                    artist1.append(text2)

        if hasattr(self, 'title'):
            for item in self.title:
                position1 = item[0]
                position2 = item[1]
                align1 = item[2]
                align2 = item[3]
                text1 = item[4]
                color1 = item[5]
                fontscale1 = item[6]
                text2 = fig1.text(
                    position1, position2, text1, horizontalalignment = align1,
                    verticalalignment = align2, transform = fig1.transFigure,
                    color = color1, fontsize = fontscale1 * self.fontsize)
            artist1.append(text2)

        gs1.tight_layout(fig1, h_pad = self.gridpad[1],
                w_pad = self.gridpad[2])

        fig1.savefig(
            filename, format = fileformat, bbox_inches = 'tight',
            bbox_extra_artists = artist1)
        plt.close()

    def __init__(
            self, data = None, kind = None, style = None, color = None,
            label = None, zorder = None, legend = None, xscale = None,
            yscale = None, xlim = None, ylim = None, xdel = None, ydel = None,
            xtick = None, ytick = None, xgrid = None, ygrid = None,
            xlabel = None, ylabel = None, xlabpos = None, ylabpos = None,
            note = None, title = None, pagesize = None, fontsize = None,
            linewidth = None, markersize = None, labelpad = None,
            tickpad = None, ticksize = None, griddim = None, gridpad = None,
            gridpos = None, frame = False):

        if data is not None:
            self.data = data
        if kind is not None:
            self.kind = kind
        if style is not None:
            self.style = style
        if color is not None:
            self.color = color
        if label is not None:
            self.label = label
        if zorder is not None:
            self.zorder = zorder
        if legend is not None:
            self.legend = legend
        if xscale is not None:
            self.xscale = xscale
        if yscale is not None:
            self.yscale = yscale
        if xlim is not None:
            self.xlim = xlim
        if ylim is not None:
            self.ylim = ylim
        if xdel is not None:
            self.xdel = xdel
        if ydel is not None:
            self.ydel = ydel
        if xtick is not None:
            self.xtick = xtick
        if ytick is not None:
            self.ytick = ytick
        if xgrid is not None:
            self.xgrid = xgrid
        if ygrid is not None:
            self.ygrid = ygrid
        if xlabel is not None:
            self.xlabel = xlabel
        if ylabel is not None:
            self.ylabel = ylabel
        if xlabpos is not None:
            self.xlabpos = xlabpos
        if ylabpos is not None:
            self.ylabpos = ylabpos
        if note is not None:
            self.note = note
        if title is not None:
            self.title = title
        if pagesize is None:
            self.pagesize = [8.0, 6.0]
        else:
            self.pagesize = pagesize
        if fontsize is None:
            self.fontsize = 18.0
        else:
            self.fontsize = fontsize
        if linewidth is None:
            self.linewidth = 1.5
        else:
            self.linewidth = linewidth
        if markersize is None:
            self.markersize = 8.0
        else:
            self.markersize = markersize
        if labelpad is None:
            self.labelpad = [6.0, 6.0]
        else:
            self.labelpad = labelpad
        if tickpad is None:
            self.tickpad = [6.0, 6.0, 6.0, 6.0]
        else:
            self.tickpad = tickpad
        if ticksize is None:
            self.ticksize = [6.0, 3.0, 6.0, 3.0]
        else:
            self.ticksize = ticksize
        if griddim is None:
            self.griddim = [1, 1]
        else:
            self.griddim = griddim
        if gridpad is None:
            self.gridpad = [0.5, 0.0, 0.0]
        else:
            self.gridpad = gridpad
        if gridpos is None:
            self.gridpos = [[0, 1, 0, 1]]
        else:
            self.gridpos = gridpos
        if frame is None:
            self.frame = False
        else:
            self.frame = frame

