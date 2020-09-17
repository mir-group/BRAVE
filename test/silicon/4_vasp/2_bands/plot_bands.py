import brave

bnd = brave.Diagram()
bnd.read('internal', ['../../../kpath/fcc.full.in'])
bnd.calc_kindex('number')
bnd.kindex[-1] -= 1
bnd.read('vasp-out', ['OUTCAR'])
bnd.calc_efermi()
bnd.set_plot('energy')
bnd.plot.ylim = [[0, 12]]
bnd.plot.pagesize = [5.0, 3.75]
bnd.plot.note = [[[0.5, 1.02, 'center', 'bottom', 'fcc Si' + 20 * ' ' + 'LDA' + 20 * ' ' + 'VASP', 'black', 1.0]]]
bnd.plot.write('matplotlib', 'png', 'silicon_bands.png')
