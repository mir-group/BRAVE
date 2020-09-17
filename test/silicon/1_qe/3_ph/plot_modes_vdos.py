import brave

bnd = brave.Diagram()
bnd.read('internal', ['../../../kpath/fcc.full.in'])
bnd.calc_kindex('density')
bnd.read('matdyn-out', ['../1_scf/1.pw.out', 'silicon.modes'])
bnd.read('matdyn-dos', ['silicon.vdos'])
bnd.set_dunit(['uc', 'cm-1'])
bnd.set_plot('energy_dos')
bnd.plot.ylim = [[0, 600], [0, 600]]
bnd.plot.pagesize = [7.0, 3.75]
bnd.plot.title = [[0.5, 1.0, 'center', 'top', 'fcc Si' + 20 * ' ' + 'LDA' + 20 * ' ' + 'QE', 'black', 1.0]]
bnd.plot.write('matplotlib', 'png', 'silicon_modes_vdos.png')
