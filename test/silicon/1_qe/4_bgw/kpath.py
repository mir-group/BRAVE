import brave
kpt = brave.Kpoint()
kpt.read('internal', ['../../../kpath/fcc.full.in'])
kpt.calc_kindex('density')
kpt.write('pw-in', ['5.pw.in.tmp'])
