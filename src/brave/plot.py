"""This module defines class Plot."""

import sys

import numpy
import matplotlib
if not 'matplotlib.pyplot' in sys.modules:
    matplotlib.use('Agg')
    import matplotlib.pyplot

import brave.common as common

class Plot(object):
    """Class for producing two-dimensional plots using
    external libraries such as matplotlib.
    """

    @property
    def data(self):
        """2D list of data curves (first dimension is plot
    number, second dimension is curve number),
    each data curve is an array of m by n floats.
    m = 2 for x, y (kind = 'plot'|'fill');
    m = 3 for x, y, s (kind = 'scatter');
    m = 4 for x, y, yerr, xerr (kind = 'errorbar');
    m = 6 for left, height, width, bottom,
    yerrdn, yerrup (kind = 'bar').
    n is number of data points (could be different
    for each data curve).
    """
        return self._data

    @data.setter
    def data(self, data):
        for group in data:
            for item in group:
                if len(item.shape) != 2 or not item.shape[0] in (2, 3, 4, 6):
                    raise ValueError(data)

        self._data = data

    @data.deleter
    def data(self):
        del self._data

    @property
    def kind(self):
        """Kinds of data curves, 2D list of strs.
    Each kind = 'plot', 'scatter', 'bar', 'errorbar',
    or 'fill'. The default is 'plot'.
    """
        return self._kind

    @kind.setter
    def kind(self, kind):
        for group in kind:
            for item in group:
                if type(item) != str:
                    raise ValueError(kind)

        self._kind = kind

    @kind.deleter
    def kind(self):
        del self._kind

    @property
    def style(self):
        """Styles of data curves, 2D list of lists.
    Each list contains styles for different elements:
    line, marker. Each linestyle = 'solid', 'dash',
    'dot', 'dashdot', (offset, on-off-dash-seq),
    or any str recognized by external library.
    Each markerstyle = 'circle', 'square', 'diamond',
    'triangle_up', 'triangle_down', 'triangle_left',
    'triangle_right', 'star', 'cross', 'plus',
    or any str recognized by external library.
    Can be set to 'None' if not used.
    The default is ['solid', 'None'].
    """
        return self._style

    @style.setter
    def style(self, style):
        for group in style:
            for subgroup in group:
                for item in subgroup:
                    if type(item) != str and type(
                            item) != tuple and item != 'None':
                        raise ValueError(style)

        self._style = style

    @style.deleter
    def style(self):
        del self._style

    @property
    def color(self):
        """Colors of data curves, 2D list of lists.
    Each list contains colors for different elements:
    line, markeredge, markerface (kind = 'plot');
    markeredge, markerface (kind = 'scatter');
    barface, baredge, errorbar (kind = 'bar');
    line, markeredge, markerface, errorbar (kind = 'errorbar');
    patchface, patchedge (kind = 'fill').
    Each color = 'black', 'red', 'green', 'blue',
    'magenta', 'yellow', 'cyan', 'white', or any
    str, float, or tuple recognized by external
    library. Can be set to 'none' if not used.
    The default is ['black', 'none', 'none'].
    """
        return self._color

    @color.setter
    def color(self, color):
        for group in color:
            for subgroup in group:
                for item in subgroup:
                    if not matplotlib.colors.is_color_like(
                            item) and item != 'none':
                        raise ValueError(color)

        self._color = color

    @color.deleter
    def color(self):
        del self._color

    @property
    def label(self):
        """Labels of data curves, 2D list of strs.
    Used for displaying legends on the plots.
    A label can be set to an empty string to
    hide it from the legend. The default is
    None.
    """
        return self._label

    @label.setter
    def label(self, label):
        for group in label:
            for item in group:
                if type(item) != str and item != None:
                    raise ValueError(label)

        self._label = label

    @label.deleter
    def label(self):
        del self._label

    @property
    def zorder(self):
        """Drawing order of data curves, 2D list of ints.
    Higher drawing order curves are drawn above
    lower drawing order curves. The default is
    the order of property data.
    """
        return self._zorder

    @zorder.setter
    def zorder(self, zorder):
        for group in zorder:
            for item in group:
                if type(item) != int and type(item) != float and item != None:
                    raise ValueError(zorder)

        self._zorder = zorder

    @zorder.deleter
    def zorder(self):
        del self._zorder

    @property
    def legend(self):
        """Legend location, list of lists.
    Each list = [loc, ncol, bbox_to_anchor] where loc
    = str is location, ncol = int is number of columns,
    bbox_to_anchor = [left, bottom, width, height] is
    position. The default is ['best', 1, None].
    """
        return self._legend

    @legend.setter
    def legend(self, legend):
        for item in legend:
            if type(item) != list and item != None:
                raise ValueError(legend)

        self._legend = legend

    @legend.deleter
    def legend(self):
        del self._legend

    @property
    def xscale(self):
        """Scaling of x-axes, list of strs. Each scaling
    = 'linear', 'log', or any string recongized
    by external library. The default is 'linear'.
    """
        return self._xscale

    @xscale.setter
    def xscale(self, xscale):
        for item in xscale:
            if type(item) != str and item != None:
                raise ValueError(xscale)

        self._xscale = xscale

    @xscale.deleter
    def xscale(self):
        del self._xscale

    @property
    def yscale(self):
        """Scaling of y-axes, list of strs. Each scaling
    = 'linear', 'log', or any string recongized
    by external library. The default is 'linear'.
    """
        return self._yscale

    @yscale.setter
    def yscale(self, yscale):
        for item in yscale:
            if type(item) != str and item != None:
                raise ValueError(yscale)

        self._yscale = yscale

    @yscale.deleter
    def yscale(self):
        del self._yscale

    @property
    def xlim(self):
        """Starting and ending values of x-axes, list of
    [xmin, xmax] where xmin and xmax are floats.
    Set to None to autoscale x-axis. The default
    is to autoscale.
    """
        return self._xlim

    @xlim.setter
    def xlim(self, xlim):
        for item in xlim:
            if type(item) != list and item != None:
                raise ValueError(xlim)

        self._xlim = xlim

    @xlim.deleter
    def xlim(self):
        del self._xlim

    @property
    def ylim(self):
        """Starting and ending values of y-axes, list of
    [ymin, ymax] where ymin and ymax are floats.
    Set to None to autoscale y-axis. The default
    is to autoscale.
    """
        return self._ylim

    @ylim.setter
    def ylim(self, ylim):
        for item in ylim:
            if type(item) != list and item != None:
                raise ValueError(ylim)

        self._ylim = ylim

    @ylim.deleter
    def ylim(self):
        del self._ylim

    @property
    def xdel(self):
        """Major and minor ticks on x-axes, list of
    [xmajor, xminor] where xmajor and xminor
    are floats. Set to None to use auto ticks
    on x-axis. The default is to use auto ticks.
    """
        return self._xdel

    @xdel.setter
    def xdel(self, xdel):
        for item in xdel:
            if type(item) != list and item != None:
                raise ValueError(xdel)

        self._xdel = xdel

    @xdel.deleter
    def xdel(self):
        del self._xdel

    @property
    def ydel(self):
        """Major and minor ticks on y-axes, list of
    [ymajor, yminor] where ymajor and yminor
    are floats. Set to None to use auto ticks
    on y-axis. The default is to use auto ticks.
    """
        return self._ydel

    @ydel.setter
    def ydel(self, ydel):
        for item in ydel:
            if type(item) != list and item != None:
                raise ValueError(ydel)

        self._ydel = ydel

    @ydel.deleter
    def ydel(self):
        del self._ydel

    @property
    def xtick(self):
        """Custom ticks on x-axes, list of lists of
    [xposition, xtext] where xposition is float
    and xtext is str. Set to None to use auto
    ticks on x-axis. The default is to use auto
    ticks.
    """
        return self._xtick

    @xtick.setter
    def xtick(self, xtick):
        for group in xtick:
            if group != None:
                for item in group:
                    if type(item) != list:
                        raise ValueError(xtick)

        self._xtick = xtick

    @xtick.deleter
    def xtick(self):
        del self._xtick

    @property
    def ytick(self):
        """Custom ticks on y-axes, list of lists of
    [yposition, ytext] where yposition is float
    and ytext is str. Set to None to use auto
    ticks on y-axis. The default is to use auto
    ticks.
    """
        return self._ytick

    @ytick.setter
    def ytick(self, ytick):
        for group in ytick:
            if group != None:
                for item in group:
                    if type(item) != list:
                        raise ValueError(ytick)

        self._ytick = ytick

    @ytick.deleter
    def ytick(self):
        del self._ytick

    @property
    def xgrid(self):
        """Grids on x-axes, list of [xgridmajorstyle,
    xgridmajorcolor, xgridmajorwidth, xgridminorstyle,
    xgridminorcolor, xgridmajorwidth], major or minor
    grids can be set to None.
    """
        return self._xgrid

    @xgrid.setter
    def xgrid(self, xgrid):
        for item in xgrid:
            if type(item) != list and item != None:
                raise ValueError(xgrid)

        self._xgrid = xgrid

    @xgrid.deleter
    def xgrid(self):
        del self._xgrid

    @property
    def ygrid(self):
        """Grids on y-axes, list of [ygridmajorstyle,
    ygridmajorcolor, ygridmajorwidth, ygridminorstyle,
    ygridminorcolor, ygridmajorwidth], major or minor
    grids can be set to None.
    """
        return self._ygrid

    @ygrid.setter
    def ygrid(self, ygrid):
        for item in ygrid:
            if type(item) != list and item != None:
                raise ValueError(ygrid)

        self._ygrid = ygrid

    @ygrid.deleter
    def ygrid(self):
        del self._ygrid

    @property
    def xlabel(self):
        """Labels of x-axes, list of strs. A label can be
    set to an empty string to hide it from the
    plot. The default is None.
    """
        return self._xlabel

    @xlabel.setter
    def xlabel(self, xlabel):
        for item in xlabel:
            if type(item) != str:
                raise ValueError(xlabel)

        self._xlabel = xlabel

    @xlabel.deleter
    def xlabel(self):
        del self._xlabel

    @property
    def ylabel(self):
        """Labels of y-axes, list of strs. A label can be
    set to an empty string to hide it from the
    plot. The default is None.
    """
        return self._ylabel

    @ylabel.setter
    def ylabel(self, ylabel):
        for item in ylabel:
            if type(item) != str:
                raise ValueError(ylabel)

        self._ylabel = ylabel

    @ylabel.deleter
    def ylabel(self):
        del self._ylabel

    @property
    def xlabpos(self):
        """Positions of x-axis labels, list of lists.
    Each list = [xposition, xalign] where
    xposition is float from 0.0 to 1.0 and
    xalign = 'center'|'right'|'left'.
    The default is None.
    """
        return self._xlabpos

    @xlabpos.setter
    def xlabpos(self, xlabpos):
        for item in xlabpos:
            if type(item) != list and item != None:
                raise ValueError(xlabpos)

        self._xlabpos = xlabpos

    @xlabpos.deleter
    def xlabpos(self):
        del self._xlabpos

    @property
    def ylabpos(self):
        """Positions of y-axis labels, list of lists.
    Each list = [yposition, yalign] where
    yposition is float from 0.0 to 1.0 and
    yalign = 'center'|'right'|'left'.
    The default is None.
    """
        return self._ylabpos

    @ylabpos.setter
    def ylabpos(self, ylabpos):
        for item in ylabpos:
            if type(item) != list and item != None:
                raise ValueError(ylabpos)

        self._ylabpos = ylabpos

    @ylabpos.deleter
    def ylabpos(self):
        del self._ylabpos

    @property
    def note(self):
        """Notes on the plots, list of lists of
    [xposition, yposition, xalign, yalign, text,
    color, fontscale] where xposition and yposition
    are floats, xalign, yalign, text, and color are
    strs, fontscale is float in units of fontsize.
    Positions are given in scaled coordinates
    which run from 0 to 1 from left to right
    and from bottom to top of the plot, xalign =
    'left'|'right' and yalign = 'bottom'|'top'.
    Set to an empty list to hide it from the plot.
    The default is None.
    """
        return self._note

    @note.setter
    def note(self, note):
        for group in note:
            for item in group:
                if type(item) != list:
                    raise ValueError(note)

        self._note = note

    @note.deleter
    def note(self):
        del self._note

    @property
    def title(self):
        """Notes on the page, list of
    [xposition, yposition, xalign, yalign, text,
    color, fontscale] where xposition and yposition
    are floats, xalign, yalign, text, and color are
    strs, fontscale is float in units of fontsize.
    Positions are given in scaled coordinates
    which run from 0 to 1 from left to right
    and from bottom to top of the page, xalign =
    'left'|'right' and yalign = 'bottom'|'top'.
    Set to an empty list to hide it from the page.
    The default is None.
    """
        return self._title

    @title.setter
    def title(self, title):
        for item in title:
            if type(item) != list:
                raise ValueError(title)

        self._title = title

    @title.deleter
    def title(self):
        del self._title

    @property
    def pagesize(self):
        """Size of the page in inches, [width, height]
    where width and height are floats.
    The default is [8.0, 6.0].
    """
        return self._pagesize

    @pagesize.setter
    def pagesize(self, pagesize):
        if type(pagesize) != list:
            raise ValueError(pagesize)

        self._pagesize = pagesize

    @pagesize.deleter
    def pagesize(self):
        del self._pagesize

    @property
    def fontsize(self):
        """Size of the font in points, float.
    The default is 18.0.
    """
        return self._fontsize

    @fontsize.setter
    def fontsize(self, fontsize):
        self._fontsize = fontsize

    @fontsize.deleter
    def fontsize(self):
        del self._fontsize

    @property
    def linewidth(self):
        """Width of the line in points, float.
    The default is 1.5.
    """
        return self._linewidth

    @linewidth.setter
    def linewidth(self, linewidth):
        self._linewidth = linewidth

    @linewidth.deleter
    def linewidth(self):
        del self._linewidth

    @property
    def markersize(self):
        """Size of the marker in points, float.
    The default is 8.0.
    """
        return self._markersize

    @markersize.setter
    def markersize(self, markersize):
        self._markersize = markersize

    @markersize.deleter
    def markersize(self):
        del self._markersize

    @property
    def labelpad(self):
        """Amount of padding for labels in points,
    [xlabelpad, ylabelpad] where xlabelpad
    and ylabelpad are floats. The default
    is [6.0, 6.0].
    """
        return self._labelpad

    @labelpad.setter
    def labelpad(self, labelpad):
        if type(labelpad) != list:
            raise ValueError(labelpad)

        self._labelpad = labelpad

    @labelpad.deleter
    def labelpad(self):
        del self._labelpad

    @property
    def tickpad(self):
        """Amount of padding for ticks in points,
    [xtickmajorpad, xtickminorpad, ytickmajorpad,
    ytickminorpad] where xtickmajorpad, xtickminorpad,
    ytickmajorpad, and ytickminorpad are floats.
    The default is [6.0, 6.0, 6.0, 6.0].
    """
        return self._tickpad

    @tickpad.setter
    def tickpad(self, tickpad):
        if type(tickpad) != list:
            raise ValueError(tickpad)

        self._tickpad = tickpad

    @tickpad.deleter
    def tickpad(self):
        del self._tickpad

    @property
    def ticksize(self):
        """Size of the ticks in points,
    [xtickmajorsize, xtickminorsize, ytickmajorsize,
    ytickminorsize] where xtickmajorsize, xtickminorsize,
    ytickmajorsize, and ytickminorsize are floats.
    The default is [6.0, 3.0, 6.0, 3.0].
    """
        return self._ticksize

    @ticksize.setter
    def ticksize(self, ticksize):
        if type(ticksize) != list:
            raise ValueError(ticksize)

        self._ticksize = ticksize

    @ticksize.deleter
    def ticksize(self):
        del self._ticksize

    @property
    def griddim(self):
        """Geometry of the grid that plots will be placed,
    [nrow, ncol] where nrow and ncol are ints.
    The default is [1, 1].
    """
        return self._griddim

    @griddim.setter
    def griddim(self, griddim):
        if type(griddim) != list:
            raise ValueError(griddim)

        self._griddim = griddim

    @griddim.deleter
    def griddim(self):
        del self._griddim

    @property
    def gridpad(self):
        """Spacing between the figure edges and the plot edges
    and height and width spacing between the plots in
    units of fontsize, [pad, h_pad, w_pad] where pad,
    h_pad, and w_pad are floats. The default is
    [0.5, 0.0, 0.0].
    """
        return self._gridpad

    @gridpad.setter
    def gridpad(self, gridpad):
        if type(gridpad) != list:
            raise ValueError(gridpad)

        self._gridpad = gridpad

    @gridpad.deleter
    def gridpad(self):
        del self._gridpad

    @property
    def gridpos(self):
        """Layout of the plots on the grid, list of
    [top, bottom, left, right] where top
    and bottom are ints in 0:nrow range,
    left and right are ints in 0:ncol
    range. The default is [[0, 1, 0, 1]].
    """
        return self._gridpos

    @gridpos.setter
    def gridpos(self, gridpos):
        for item in gridpos:
            if type(item) != list:
                raise ValueError(gridpos)

        self._gridpos = gridpos

    @gridpos.deleter
    def gridpos(self):
        del self._gridpos

    @property
    def frame(self):
        """Display legend frame, bool.
    The default is False.
    """
        return self._frame

    @frame.setter
    def frame(self, frame):
        if type(frame) != bool:
            raise ValueError(frame)

        self._frame = frame

    @frame.deleter
    def frame(self):
        del self._frame

    def write(self, libraryname, fileformat, filename):
        """Method for writing plot to file.

    libraryname     fileformat           filename
    -----------     ----------           --------
    'matplotlib'    'png'|'eps'|'pdf'    'plot.ext'
        """

        if libraryname.lower() == 'matplotlib':
            self._writep_matplotlib(fileformat, filename)
        else:
            raise ValueError(libraryname)

    def __init__(
            self, data=None, kind=None, style=None, color=None,
            label=None, zorder=None, legend=None, xscale=None, yscale=None,
            xlim=None, ylim=None, xdel=None, ydel=None, xtick=None,
            ytick=None, xgrid=None, ygrid=None, xlabel=None,
            ylabel=None, xlabpos=None, ylabpos=None,
            note=None, title=None, pagesize=[8.0, 6.0],
            fontsize=18.0, linewidth=1.5, markersize=8.0,
            labelpad=[6.0, 6.0], tickpad=[6.0, 6.0, 6.0, 6.0],
            ticksize=[6.0, 3.0, 6.0, 3.0], griddim=[1, 1],
            gridpad=[0.5, 0.0, 0.0], gridpos=[[0, 1, 0, 1]],
            frame=False):

        if data != None:
            self.data = data
        if kind != None:
            self.kind = kind
        if style != None:
            self.style = style
        if color != None:
            self.color = color
        if label != None:
            self.label = label
        if zorder != None:
            self.zorder = zorder
        if legend != None:
            self.legend = legend
        if xscale != None:
            self.xscale = xscale
        if yscale != None:
            self.yscale = yscale
        if xlim != None:
            self.xlim = xlim
        if ylim != None:
            self.ylim = ylim
        if xdel != None:
            self.xdel = xdel
        if ydel != None:
            self.ydel = ydel
        if xtick != None:
            self.xtick = xtick
        if ytick != None:
            self.ytick = ytick
        if xgrid != None:
            self.xgrid = xgrid
        if ygrid != None:
            self.ygrid = ygrid
        if xlabel != None:
            self.xlabel = xlabel
        if ylabel != None:
            self.ylabel = ylabel
        if xlabpos != None:
            self.xlabpos = xlabpos
        if ylabpos != None:
            self.ylabpos = ylabpos
        if note != None:
            self.note = note
        if title != None:
            self.title = title
        self.pagesize = pagesize
        self.fontsize = fontsize
        self.linewidth = linewidth
        self.markersize = markersize
        self.labelpad = labelpad
        self.tickpad = tickpad
        self.ticksize = ticksize
        self.griddim = griddim
        self.gridpad = gridpad
        self.gridpos = gridpos
        self.frame = frame

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
                    if _label != '' and _label != None:
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
                    if type(_aa) == str:
                        if _ii == 0:
                            if _aa.lower() in _linestyle_dict:
                                _bb = _linestyle_dict[_aa.lower()]
                        else:
                            if _aa.lower() in _markerstyle_dict:
                                _bb = _markerstyle_dict[_aa.lower()]
                    _style_matplotlib.append(_bb)

                _color_matplotlib = []
                for _cc in _color:
                    _dd = _cc
                    if type(_cc) == str:
                        if _cc.lower() in _color_dict:
                            _dd = _color_dict[_cc.lower()]
                    _color_matplotlib.append(_dd)

                if _kind.lower() == 'plot':
                    _curve = _plot.plot(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            linestyle = _style_matplotlib[0],
                            marker = _style_matplotlib[1],
                            color = _color_matplotlib[0],
                            markeredgecolor = _color_matplotlib[1],
                            markerfacecolor = _color_matplotlib[2],
                            label = _label, zorder = _zorder)

                elif _kind.lower() == 'scatter':
                    _curve = _plot.scatter(
                            self.data[iplot][idata][0],
                            self.data[iplot][idata][1],
                            s = self.data[iplot][idata][2],
                            marker = _style_matplotlib[1],
                            edgecolors = _color_matplotlib[0],
                            color = _color_matplotlib[1],
                            label = _label, zorder = _zorder)

                elif _kind.lower() == 'bar':
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

                elif _kind.lower() == 'errorbar':
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

                elif _kind.lower() == 'fill':
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
                    if self.legend[iplot] != None:
                        if type(self.legend[iplot][0]) == str:
                            _loc = self.legend[iplot][0]
                        else:
                            _loc = 'best'
                        if type(self.legend[iplot][1]) == int:
                            _ncol = self.legend[iplot][1]
                        else:
                            _ncol = 1
                        if type(self.legend[iplot][2]) == list:
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
                if self.xlim[iplot] != None:
                    _plot.set_xlim(self.xlim[iplot][0], self.xlim[iplot][1])

            if hasattr(self, 'ylim'):
                if self.ylim[iplot] != None:
                    _plot.set_ylim(self.ylim[iplot][0], self.ylim[iplot][1])

            _xmin, _xmax, _ymin, _ymax = _plot.axis()

            if hasattr(self, 'xdel'):
                if self.xdel[iplot] != None:
                    _plot.xaxis.set_major_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.xdel[iplot][0]))
                    _plot.xaxis.set_minor_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.xdel[iplot][1]))

            if hasattr(self, 'ydel'):
                if self.ydel[iplot] != None:
                    _plot.yaxis.set_major_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.ydel[iplot][0]))
                    _plot.yaxis.set_minor_locator(
                            matplotlib.ticker.MultipleLocator(
                            self.ydel[iplot][1]))

            if hasattr(self, 'xtick'):
                if self.xtick[iplot] != None:
                    _xtickposition = []
                    _xticklabel = []
                    for item in self.xtick[iplot]:
                        _xtickposition.append(item[0])
                        _xticklabel.append(item[1])
                    _plot.set_xticks(_xtickposition)
                    _plot.set_xticklabels(_xticklabel)

            if hasattr(self, 'ytick'):
                if self.ytick[iplot] != None:
                    _ytickposition = []
                    _yticklabel = []
                    for item in self.ytick[iplot]:
                        _ytickposition.append(item[0])
                        _yticklabel.append(item[1])
                    _plot.set_yticks(_ytickposition)
                    _plot.set_yticklabels(_yticklabel)

            if hasattr(self, 'xgrid'):
                if self.xgrid[iplot] != None:
                    _style = self.xgrid[iplot][0]
                    _color = self.xgrid[iplot][1]
                    _width = self.xgrid[iplot][2]
                    if _style != None and _color != None and _width != None:
                        _style_matplotlib = _style
                        if type(_style) == str:
                            if _style.lower() in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[
                                        _style.lower()]
                        _color_matplotlib = _color
                        if type(_color) == str:
                            if _color.lower() in _color_dict:
                                _color_matplotlib = _color_dict[_color.lower()]
                        _plot.grid(
                                b = True, which = 'major', axis = 'x',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

                    _style = self.xgrid[iplot][3]
                    _color = self.xgrid[iplot][4]
                    _width = self.xgrid[iplot][5]
                    if _style != None and _color != None and _width != None:
                        _style_matplotlib = _style
                        if type(_style) == str:
                            if _style.lower() in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[
                                        _style.lower()]
                        _color_matplotlib = _color
                        if type(_color) == str:
                            if _color.lower() in _color_dict:
                                _color_matplotlib = _color_dict[_color.lower()]
                        _plot.grid(
                                b = True, which = 'minor', axis = 'x',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

            if hasattr(self, 'ygrid'):
                if self.ygrid[iplot] != None:
                    _style = self.ygrid[iplot][0]
                    _color = self.ygrid[iplot][1]
                    _width = self.ygrid[iplot][2]
                    if _style != None and _color != None and _width != None:
                        _style_matplotlib = _style
                        if type(_style) == str:
                            if _style.lower() in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[
                                        _style.lower()]
                        _color_matplotlib = _color
                        if type(_color) == str:
                            if _color.lower() in _color_dict:
                                _color_matplotlib = _color_dict[_color.lower()]
                        _plot.grid(
                                b = True, which = 'major', axis = 'y',
                                linestyle = _style_matplotlib,
                                color = _color_matplotlib,
                                linewidth = _width)

                    _style = self.ygrid[iplot][3]
                    _color = self.ygrid[iplot][4]
                    _width = self.ygrid[iplot][5]
                    if _style != None and _color != None and _width != None:
                        _style_matplotlib = _style
                        if type(_style) == str:
                            if _style.lower() in _linestyle_dict:
                                _style_matplotlib = _linestyle_dict[
                                        _style.lower()]
                        _color_matplotlib = _color
                        if type(_color) == str:
                            if _color.lower() in _color_dict:
                                _color_matplotlib = _color_dict[_color.lower()]
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
                        if self.xlabpos[iplot] != None:
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
                        if self.ylabpos[iplot] != None:
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

