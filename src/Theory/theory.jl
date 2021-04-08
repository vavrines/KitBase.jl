# ============================================================
# Theories
# ============================================================

export prim_conserve, 
       conserve_prim, 
       mixture_prim_conserve, 
       mixture_conserve_prim, 
       em_coefficients, 
       advection_flux, 
       burgers_flux, 
       euler_flux, 
       euler_jacobi
export gauss_moments, 
       mixture_gauss_moments,
       moments_conserve,
       diatomic_moments_conserve,
       mixture_moments_conserve,
       pdf_slope,
       mixture_pdf_slope,
       moments_conserve_slope,
       mixture_moments_conserve_slope,
       discrete_moments,
       pressure,
       stress,
       heat_flux,
       maxwellian,
       maxwellian!,
       mixture_maxwellian,
       mixture_maxwellian!,
       shakhov,
       shakhov!,
       rykov!,
       reduce_distribution,
       full_distribution,
       ref_vhs_vis,
       vhs_collision_time,
       rykov_zr,
       aap_hs_collision_time,
       aap_hs_prim,
       aap_hs_diffeq!,
       shift_pdf!,
       hs_boltz_kn,
       kernel_mode,
       boltzmann_fft,
       boltzmann_fft!,
       boltzmann_ode!,
       bgk_ode!,
       chapman_enskog
export heat_capacity_ratio,
       sound_speed
export sample_maxwell,
       next_collision_time,
       boundary_time

include("theory_continuum.jl")
include("theory_plasma.jl")
include("theory_maxwellian.jl")
include("theory_shakhov.jl")
include("theory_diatomic.jl")
include("theory_atom.jl")
include("theory_pdf.jl")
include("theory_fsm.jl")
include("theory_moments_pure.jl")
include("theory_moments_mixture.jl")
include("theory_thermo.jl")
include("theory_particle.jl")