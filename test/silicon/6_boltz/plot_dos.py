import brave

bnd = brave.Diagram()
bnd.read('boltztrap-dos', ['silicon.intrans', 'silicon.transdos'])
bnd.set_dunit(['uc', 'ev'])
bnd.set_plot('dos')
bnd.plot.ylim = [[-6, 6]]
bnd.plot.pagesize = [2.5, 3.75]
bnd.plot.note = [[[0.5, 1.02, 'center', 'bottom', 'fcc Si' + 8 * ' ' + 'LDA' + 8 * ' ' + 'QE', 'black', 1.0]]]
bnd.plot.write('matplotlib', 'png', 'silicon_dos.png')
