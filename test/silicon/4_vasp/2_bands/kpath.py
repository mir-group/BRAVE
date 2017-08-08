import brave
kpt = brave.Kpoint()
kpt.read('internal', ['../../../kpath/fcc.full.in'])
kpt.calc_kindex('number')
kpt.write('vasp-kpt', ['KPOINTS'])
