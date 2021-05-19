# ============================================================
# Abstract Types
# ============================================================

abstract type AbstractPhysicalSpace end
abstract type AbstractStructPhysicalSpace <: AbstractPhysicalSpace end
abstract type AbstractUnstructPhysicalSpace <: AbstractPhysicalSpace end
abstract type AbstractPhysicalSpace1D <: AbstractPhysicalSpace end
abstract type AbstractPhysicalSpace2D <: AbstractPhysicalSpace end
abstract type AbstractPhysicalSpace3D <: AbstractPhysicalSpace end
abstract type AbstractVelocitySpace end
abstract type AbstractVelocitySpace1D <: AbstractVelocitySpace end
abstract type AbstractVelocitySpace2D <: AbstractVelocitySpace end
abstract type AbstractVelocitySpace3D <: AbstractVelocitySpace end
abstract type AbstractSetup end
abstract type AbstractProperty end
abstract type AbstractCondition end
abstract type AbstractSolverSet end

abstract type AbstractControlVolume end
abstract type AbstractInterface end
abstract type AbstractControlVolume1D <: AbstractControlVolume end
abstract type AbstractControlVolume2D <: AbstractControlVolume end
abstract type AbstractUnstructControlVolume <: AbstractControlVolume end
abstract type AbstractInterface1D <: AbstractInterface end
abstract type AbstractInterface2D <: AbstractInterface end

abstract type AbstractSolution end
abstract type AbstractFlux end
abstract type AbstractSolution1D <: AbstractSolution end
abstract type AbstractSolution2D <: AbstractSolution end
abstract type AbstractFlux1D <: AbstractFlux end
abstract type AbstractFlux2D <: AbstractFlux end

abstract type AbstractParticle end
abstract type AbstractParticle1D <: AbstractParticle end
abstract type AbstractParticle2D <: AbstractParticle end
