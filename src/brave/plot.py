"""This module defines class Plot."""

import sys

import numpy
import matplotlib
if not 'matplotlib.pyplot' in sys.modules:
    matplotlib.use('Agg')
    import matplotlib.pyplot

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
                item, numpy.ndarray) and item.dtype == numpy.dtype(
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
        if not all(isinstance(item, str) for group in value for item in group):
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
        """Method for writing plot to file.

    libraryname     fileformat           filename
    -----------     ----------           --------
    'matplotlib'    'png'|'eps'|'pdf'    'plot.ext'
        """

        if libraryname == 'matplotlib':
            self._writep_matplotlib(fileformat, filename)
        else:
            raise ValueError(libraryname)

    def _writep_matplotlib(self, fileformat, filename):

        _linestyle_dict = {
                'solid': '-',
                'dash': '--',
                'dot': ':',
                'dashdot': '-.'}

        _markerstyle_dict = {
                'circle': 'o',
                'square': 's',
                'diamond': 'D',
                'triangle_up': '^',
                'triangle_down': 'v',
                'triangle_left': '<',
                'triangle_right': '>',
                'star': '*',
                'cross': 'x',
                'plus': '+'}

        _color_dict = {
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

        _figure = matplotlib.pyplot.figure(figsize = (
                self.pagesize[0], self.pagesize[1]))
        _grid = matplotlib.gridspec.GridSpec(self.griddim[0], self.griddim[1])

        _artist_list = []
        for iplot in range(len(self.data)):
            _plot = _figure.add_subplot(_grid[
                    self.gridpos[iplot][0]:self.gridpos[iplot][1],
                    self.gridpos[iplot][2]:self.gridpos[iplot][3]],
                    axisbelow = True)

            _is_legend = False
            for idata in range(len(self.data[iplot])):
                if hasattr(self, 'kind'):
                    _kind = self.kind[iplot][idata]
                else:
                    _kind = 'plot'

                if hasattr(self, 'style'):
                    _style = self.style[iplot][idata]
                else:
                    _style = ['solid', 'None']

                if hasattr(self, 'color'):
                    _color = self.color[iplot][idata]
                else:
                    _color = ['black', 'none', 'none']

                if hasattr(self, 'label'):
                    _label = self.label[iplot][idata]
                    if _label != '' and _label is not None:
                        _is_legend = True
                else:
                    _label = ''

                if hasattr(self, 'zorder'):
                    _zorder = self.zorder[iplot][idata]
                else:
                    _zorder = 2

                _style_matplotlib = []
                for _ii in range(2):
                    _aa = _style[_ii]
                    _bb = _aa
                    if type(_aa) is str:
                        if _ii == 0:
                            if _aa in _linestyle_dict:
                                _bb = _linestyle_dict[_aa]
                        else:
                            if _aa in _markerstyle_dict:
                                _bb = _markerstyle_dict[_aa]
                    _style_matplotlib.append(_bb)

                _color_matplotlib = []
                for _cc in _color:
                    _dd = _cc
                    if type(_cc) is str:
                        if _cc in _color_dict:
                            _dd = _color_dict[_cc]
                    _color_matplotlib.append(_dd)

                if _kind == 'plot':
                    _curve = _plot.plot(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            linestyle = _style_matplotlib[0],
                            marker = _style_matplotlib[1],
                            color = _color_matplotlib[0],
                            markeredgecolor = _color_matplotlib[1],
                            markerfacecolor = _color_matplotlib[2],
                            label = _label, zorder = _zorder)

                elif _kind == 'scatter':
                    _curve = _plot.scatter(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            s = self.data[iplot][idata][2],
                            marker = _style_matplotlib[1],
                            edgecolors = _color_matplotlib[0],
                            color = _color_matplotlib[1],
                            label = _label, zorder = _zorder)

                elif _kind == 'bar':
                    if max(
                            numpy.amax(self.data[iplot][idata][4]),
                            numpy.amax(self.data[iplot][idata][5])
                            ) > common.EPS12:
                        yerr = [
                                self.data[iplot][idata][4],
                                self.data[iplot][idata][5]]
                    else:
                        yerr = None
                    _curve = _plot.bar(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            width = self.data[iplot][idata][2],
                            bottom = self.data[iplot][idata][3],
                            yerr = yerr,
                            color = _color_matplotlib[0],
                            edgecolor = _color_matplotlib[1],
                            ecolor = _color_matplotlib[2],
                            linewidth = self.linewidth,
                            label = _label, zorder = _zorder)

                elif _kind == 'errorbar':
                    if numpy.amax(self.data[iplot][idata][2]) > common.EPS12:
                        yerr = self.data[iplot][idata][2]
                    else:
                        yerr = None
                    if numpy.amax(self.data[iplot][idata][3]) > common.EPS12:
                        xerr = self.data[iplot][idata][3]
                    else:
                        xerr = None
                    _curve = _plot.errorbar(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            yerr = yerr, xerr = xerr,
                            linestyle = _style_matplotlib[0],
                            marker = _style_matplotlib[1],
                            color = _color_matplotlib[0],
                            markeredgecolor = _color_matplotlib[1],
                            markerfacecolor = _color_matplotlib[2],
                            ecolor = _color_matplotlib[3],
                            label = _label, zorder = _zorder)

                elif _kind == 'fill':
                    _curve = _plot.fill(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            linestyle = _style_matplotlib[0],
                            facecolor = _color_matplotlib[0],
                            edgecolor = _color_matplotlib[1],
                            label = _label, zorder = _zorder)

                else:
                    raise ValueError(_kind)

            if _is_legend:
                if hasattr(self, 'legend'):
                    if self.legend[iplot] is not None:
                        if type(self.legend[iplot][0]) is str:
                            _loc = self.legend[iplot][0]
                        else:
                            _loc = 'best'
                        if type(self.legend[iplot][1]) is int:
                            _ncol = self.legend[iplot][1]
                        else:
                            _ncol = 1
                        if type(self.legend[iplot][2]) is list:
                            _borderaxespad = 0.0
                            _bbox_to_anchor = self.legend[iplot][2]
                            _bbox_transform = _figure.transFigure
                        else:
                            _borderaxespad = 0.5
                            _bbox_to_anchor = None
                            _bbox_transform = _plot.transAxes
                    else:
                        _loc = 'best'
                        _ncol = 1
                        _borderaxespad = 0.5
                        _bbox_to_anchor = None
                        _bbox_transform = _plot.transAxes
                else:
                    _loc = 'best'
                    _ncol = 1
                    _borderaxespad = 0.5
                    _bbox_to_anchor = None
                    _bbox_transform = _plot.transAxes
                _legend = _plot.legend(
                        loc = _loc, ncol = _ncol,
                        borderaxespad = _borderaxespad,
                        bbox_to_anchor = _bbox_to_anchor,
                        bbox_transform = _bbox_transform)
                if self.frame:
                    _legend.get_frame().set_linewidth(self.linewidth)

            if hasattr(self, 'xscale'):
                _plot.set_xscale(self.xscale[iplot])

            if hasattr(self, 'yscale'):
                _plot.set_yscale(self.yscale[iplot])

            _plot.autoscale(enable = True, axis = 'both', tight = False)

            if hasattr(self, 'xlim'):
                if self.xlim[iplot] is not None:
                    _plot.set_xlim(self.xlim[iplot][0], self.xlim[iplot][1])

            if hasattr(self, 'ylim'):
                if self.ylim[iplot] is not None:
                    _plot.set_ylim(self.ylim[iplot][0], self.ylim[iplot][1])

            _xmin, _xmax, _ymin, _ymax = _plot.axis()

            if hasattr(self, 'xdel'):
                if self.xdel[iplot] is not None:
                    _plot.xaxis.set_major_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.xdel[iplot][0]))
                    _plot.xaxis.set_minor_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.xdel[iplot][1]))

            if hasattr(self, 'ydel'):
                if self.ydel[iplot] is not None:
                    _plot.yaxis.set_major_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.ydel[iplot][0]))
                    _plot.yaxis.set_minor_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.ydel[iplot][1]))

            if hasattr(self, 'xtick'):
                if self.xtick[iplot] is not None:
                    _xtickposition = []
                    _xticklabel = []
                    for item in self.xtick[iplot]:
                        _xtickposition.append(item[0])
                        _xticklabel.append(item[1])
                    _plot.set_xticks(_xtickposition)
                    _plot.set_xticklabels(_xticklabel)

            if hasattr(self, 'ytick'):
                if self.ytick[iplot] is not None:
                    _ytickposition = []
                    _yticklabel = []
                    for item in self.ytick[iplot]:
                        _ytickposition.append(item[0])
                        _yticklabel.append(item[1])
                    _plot.set_yticks(_ytickposition)
                    _plot.set_yticklabels(_yticklabel)

            if hasattr(self, 'xgrid'):
                if self.xgrid[iplot][0] is not None:
                    _style = self.xgrid[iplot][0][0]
                    _color = self.xgrid[iplot][0][1]
                    _width = self.xgrid[iplot][0][2]
                    if not None in (_style, _color, _width):
                        _style_matplotlib = _style
                        if type(_style) is str:
                            if _style in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[_style]
                        _color_matplotlib = _color
                        if type(_color) is str:
                            if _color in _color_dict:
                                _color_matplotlib = _color_dict[_color]
                        _plot.grid(
                                b = True, which = 'major', axis = 'x',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

                if self.xgrid[iplot][1] is not None:
                    _style = self.xgrid[iplot][1][0]
                    _color = self.xgrid[iplot][1][1]
                    _width = self.xgrid[iplot][1][2]
                    if not None in (_style, _color, _width):
                        _style_matplotlib = _style
                        if type(_style) is str:
                            if _style in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[_style]
                        _color_matplotlib = _color
                        if type(_color) is str:
                            if _color in _color_dict:
                                _color_matplotlib = _color_dict[_color]
                        _plot.grid(
                                b = True, which = 'minor', axis = 'x',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

            if hasattr(self, 'ygrid'):
                if self.ygrid[iplot][0] is not None:
                    _style = self.ygrid[iplot][0][0]
                    _color = self.ygrid[iplot][0][1]
                    _width = self.ygrid[iplot][0][2]
                    if not None in (_style, _color, _width):
                        _style_matplotlib = _style
                        if type(_style) is str:
                            if _style in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[_style]
                        _color_matplotlib = _color
                        if type(_color) is str:
                            if _color in _color_dict:
                                _color_matplotlib = _color_dict[_color]
                        _plot.grid(
                                b = True, which = 'major', axis = 'y',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

                if self.ygrid[iplot][1] is not None:
                    _style = self.ygrid[iplot][1][0]
                    _color = self.ygrid[iplot][1][1]
                    _width = self.ygrid[iplot][1][2]
                    if not None in (_style, _color, _width):
                        _style_matplotlib = _style
                        if type(_style) is str:
                            if _style in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[_style]
                        _color_matplotlib = _color
                        if type(_color) is str:
                            if _color in _color_dict:
                                _color_matplotlib = _color_dict[_color]
                        _plot.grid(
                                b = True, which = 'minor', axis = 'y',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

            if hasattr(self, 'xlabel'):
                if self.xlabel[iplot] == '':
                    _xticklabels = []
                    for _xticklabel in _plot.get_xticklabels():
                        _xticklabels.append('')
                    _plot.set_xticklabels(_xticklabels)
                else:
                    _xposition = 0.5
                    _xalign = 'center'
                    if hasattr(self, 'xlabpos'):
                        if self.xlabpos[iplot] is not None:
                            _xposition = self.xlabpos[iplot][0]
                            _xalign = self.xlabpos[iplot][1]
                    _xlabel = _plot.set_xlabel(
                            self.xlabel[iplot], labelpad = self.labelpad[0],
                            multialignment = 'center', x = _xposition,
                            ha = _xalign)
                    _artist_list.append(_xlabel)

            if hasattr(self, 'ylabel'):
                if self.ylabel[iplot] == '':
                    _yticklabels = []
                    for _yticklabel in _plot.get_yticklabels():
                        _yticklabels.append('')
                    _plot.set_yticklabels(_yticklabels)
                else:
                    _yposition = 0.5
                    _yalign = 'center'
                    if hasattr(self, 'ylabpos'):
                        if self.ylabpos[iplot] is not None:
                            _yposition = self.ylabpos[iplot][0]
                            _yalign = self.ylabpos[iplot][1]
                    _ylabel = _plot.set_ylabel(
                            self.ylabel[iplot], labelpad = self.labelpad[1],
                            multialignment = 'center', y = _yposition,
                            ha = _yalign)
                    _artist_list.append(_ylabel)

            if hasattr(self, 'note'):
                for item in self.note[iplot]:
                    _xposition = item[0]
                    _yposition = item[1]
                    _xalign = item[2]
                    _yalign = item[3]
                    _text = item[4]
                    _color = item[5]
                    _fontscale = item[6]
                    _note = _plot.text(
                            _xposition, _yposition, _text,
                            horizontalalignment = _xalign,
                            verticalalignment = _yalign,
                            transform = _plot.transAxes,
                            color = _color,
                            fontsize = _fontscale * self.fontsize)
                    _artist_list.append(_note)

        if hasattr(self, 'title'):
            for item in self.title:
                _xposition = item[0]
                _yposition = item[1]
                _xalign = item[2]
                _yalign = item[3]
                _text = item[4]
                _color = item[5]
                _fontscale = item[6]
                _title = _figure.text(
                        _xposition, _yposition, _text,
                        horizontalalignment = _xalign,
                        verticalalignment = _yalign,
                        transform = _figure.transFigure,
                        color = _color,
                        fontsize = _fontscale * self.fontsize)
            _artist_list.append(_title)

        _grid.tight_layout(_figure, h_pad = self.gridpad[1],
                w_pad = self.gridpad[2])

        _figure.savefig(
                filename, format = fileformat, bbox_inches = 'tight',
                bbox_extra_artists = _artist_list)
        matplotlib.pyplot.close()

    def __init__(
            self, data=None, kind=None, style=None, color=None, label=None,
            zorder=None, legend=None, xscale=None, yscale=None, xlim=None,
            ylim=None, xdel=None, ydel=None, xtick=None, ytick=None,
            xgrid=None, ygrid=None, xlabel=None, ylabel=None, xlabpos=None,
            ylabpos=None, note=None, title=None, pagesize=None, fontsize=None,
            linewidth=None, markersize=None, labelpad=None, tickpad=None,
            ticksize=None, griddim=None, gridpad=None, gridpos=None,
            frame=False):

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

