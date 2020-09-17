import brave

bnd = brave.Diagram()
bnd.read('internal', ['../../kpath/fcc.full.in'])
bnd.calc_kindex('density')
bnd.calc_kpoint()
bnd.read('wannier-out', ['silicon.win', 'silicon-bands.dat'])
bnd.nelec = 8.0
bnd.calc_efermi()
bnd.set_plot('energy')
bnd.plot.ylim = [[0, 12]]
bnd.plot.pagesize = [5.0, 3.75]
bnd.plot.note = [[[0.5, 1.02, 'center', 'bottom', 'fcc Si' + 20 * ' ' + 'LDA' + 20 * ' ' + 'Wannier90', 'black', 1.0]]]
bnd.plot.write('matplotlib', 'png', 'silicon_bands.png')
