using Test, KitBase

cd(@__DIR__)
include("test_data.jl")
include("test_struct.jl")
include("test_io.jl")
include("test_math.jl")
include("test_geo.jl")
include("test_theory.jl")
include("test_moments.jl")
include("test_quadrature.jl")
include("test_config.jl")
include("test_flux.jl")
include("test_reconstruction.jl")
include("test_solver_scalar.jl")
include("test_solver_1d.jl")
include("test_solver_2d.jl")
include("test_step.jl")
include("test_particle.jl")
include("test_unstructure.jl")

ks = SolverSet("config_1d1f1v.txt")
ctr, face = init_fvm(ks, ks.ps)

ks.ib.fw[ks.ps.x[1]]



ks.ib.fw(0.0)
ks.ib.fw(1.0)

ks.ib.bcL

ks.ib.bcR