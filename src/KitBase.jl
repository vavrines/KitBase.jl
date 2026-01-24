"""
KitBase.jl: The lightweight module of solution algorithms in Kinetic.jl

Copyright (c) 2020-2025 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitBase

if VERSION < v"1.3"
    @warn "To use all the features of Kinetic, please upgrade to Julia 1.3 or newer."
end

const KB = KitBase

import BSON
import CSV
import JLD2
import NonlinearSolve
import SciMLNLSolve: NLSolveJL

using Base.Threads: @threads
using CUDA
using Dates
using Distributions
using FastGaussQuadrature
using FFTW
using FileIO
using ForwardDiff
using LinearAlgebra
using MultivariatePolynomials
using OffsetArrays
using Optim
using Parameters
using Printf: @printf
using RecipesBase
using Reexport
using SpecialFunctions
using StaticArrays
using StructArrays
using TypedPolynomials
using WriteVTK

@reexport using FiniteMesh
using FiniteMesh.DocStringExtensions
using FiniteMesh.ProgressMeter

include("Data/data.jl")
include("Macro/macro.jl")
include("Struct/struct.jl")
include("IO/io.jl")
include("Math/math.jl")
include("Geometry/geometry.jl")
include("Theory/theory.jl")
include("Phase/phase.jl")
include("Reconstruction/reconstruction.jl")
include("Flux/flux.jl")
include("AD/ad.jl")
include("Config/config.jl")
include("Boundary/boundary.jl")
include("Solver/solver.jl")

export KB

end
