import brave
kpt = brave.Kpoint()
kpt.read('internal', ['../../../kpath/fcc.full.in'])
kpt.calc_kindex('number')
kpt.write('lapw-kpt', ['silicon.klist_band'], lapwkunit='cartesian')
