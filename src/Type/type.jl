# ============================================================
# Type Hierarchy
# ============================================================

export AbstractSolverSet
export AbstractSetup, AbstractProperty, AbstractCondition
export AbstractStructPhysicalSpace, AbstractUnstructPhysicalSpace
export AbstractPhysicalSpace, AbstractPhysicalSpace1D, AbstractPhysicalSpace2D
export AbstractVelocitySpace, AbstractVelocitySpace1D, AbstractVelocitySpace2D
export AbstractControlVolume, AbstractControlVolume1D, AbstractControlVolume2D
export AbstractUnstructControlVolume
export AbstractInterface, AbstractInterface1D, AbstractInterface2D
export AbstractSolution, AbstractSolution1D, AbstractSolution2D
export AbstractFlux, AbstractFlux1D, AbstractFlux2D
export AbstractParticle, AbstractParticle1D, AbstractParticle2D

export Setup
export Scalar, Radiation, Gas, DiatomicGas, Mixture, Plasma1D, Plasma2D
export IB, IB1F, IB2F, IB3F, IB4F
export ControlVolume1D, ControlVolume1D1F, ControlVolume1D2F, ControlVolume1D3F, ControlVolume1D4F
export ControlVolume2D, ControlVolume2D1F, ControlVolume2D2F, ControlVolume2D3F
export ControlVolumeUS, ControlVolumeUS1F, ControlVolumeUS2F
export Interface1D, Interface1D1F, Interface1D2F, Interface1D3F, Interface1D4F
export Interface2D, Interface2D1F, Interface2D2F
export Solution1D, Solution1D1F, Solution1D2F
export Solution2D, Solution2D1F, Solution2D2F
export Flux1D, Flux1D1F, Flux1D2F
export Flux2D, Flux2D1F, Flux2D2F
export Particle, Particle1D, Particle2D
export ControlVolumeParticle1D, ControlVolumeParticle2D

include("abstract.jl")
include("struct_general.jl")
include("struct_ib.jl")
include("struct_ctr.jl")
include("struct_face.jl")
include("struct_sol.jl")
include("struct_flux.jl")
include("struct_ptc.jl")
