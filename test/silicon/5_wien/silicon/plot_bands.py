import brave

bnd = brave.Diagram()
bnd.read('internal', ['../../../kpath/fcc.full.in'])
bnd.calc_kindex('number')
bnd.read('lapw-out', ['silicon.output1'])
bnd.set_eunit('ev')
bnd.calc_efermi()
bnd.set_plot('energy')
bnd.plot.ylim = [[0, 12]]
bnd.plot.pagesize = [5.0, 3.75]
bnd.plot.note = [[[0.5, 1.02, 'center', 'bottom', 'fcc Si' + 20 * ' ' + 'LDA' + 20 * ' ' + 'WIEN2k', 'black', 1.0]]]
bnd.plot.write('matplotlib', 'png', 'silicon_bands.png')
