import brave

bnd = brave.Diagram()
bnd.read('matdyn-dos', ['silicon.vdos'])
bnd.set_plot('dos')
bnd.plot.ylim = [[0, 600]]
bnd.plot.pagesize = [2.5, 3.75]
bnd.plot.note = [[[0.5, 1.02, 'center', 'bottom', 'fcc Si' + 8 * ' ' + 'LDA' + 8 * ' ' + 'QE', 'black', 1.0]]]
bnd.plot.write('matplotlib', 'png', 'silicon_vdos.png')
