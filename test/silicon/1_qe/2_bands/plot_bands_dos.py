import brave

bnd = brave.Diagram()
bnd.read('internal', ['../../../kpath/fcc.full.in'])
bnd.calc_kindex('density')
bnd.read('pw-out', ['1.pw.out'])
bnd.calc_efermi()
bnd.read('boltztrap-dos', ['../../6_boltz/silicon.intrans', '../../6_boltz/silicon.transdos'])
bnd.set_dunit(['uc', 'ev'])
bnd.dos[0, :] += bnd.efermi
bnd.set_plot('energy_dos')
bnd.plot.ylim = [[0, 12], [0, 12]]
bnd.plot.pagesize = [7.0, 3.75]
bnd.plot.title = [[0.5, 1.0, 'center', 'top', 'fcc Si' + 20 * ' ' + 'LDA' + 20 * ' ' + 'QE', 'black', 1.0]]
bnd.plot.write('matplotlib', 'png', 'silicon_bands_dos.png')
