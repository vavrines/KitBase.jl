using KitBase.OffsetArrays

KitBase.static_array(zeros(0:2))
KitBase.static_array(zeros(0:2, 0:1))
KitBase.static_array(zeros(0:2, 0:1, 0:1))
KitBase.static_array(zeros(0:2, 0:1, 0:1, 0:1))

KitBase.collect_run(`ls`)
KitBase.generate_vars(Dict(:a => 1))
