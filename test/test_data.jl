using KitBase.OffsetArrays

KB.static_array(zeros(0:2))
KB.static_array(zeros(0:2, 0:1))
KB.static_array(zeros(0:2, 0:1, 0:1))
KB.static_array(zeros(0:2, 0:1, 0:1, 0:1))

KB.collect_run(`ls`)
KB.generate_vars(Dict(:a => 1))

KB.extract_last(rand(2, 2), 1; mode = :view)
KB.extract_last(rand(2, 2, 2), 1; mode = :view)
KB.extract_last(rand(2, 2, 2, 2), 1; mode = :view)
KB.extract_last(rand(2, 2, 2, 2, 2), 1; mode = :view)
KB.extract_last(rand(2, 2), 1; mode = :copy)
KB.extract_last(rand(2, 2, 2), 1; mode = :copy)
KB.extract_last(rand(2, 2, 2, 2), 1; mode = :copy)
KB.extract_last(rand(2, 2, 2, 2, 2), 1; mode = :copy)
