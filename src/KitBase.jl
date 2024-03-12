"""
KitBase.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2023 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitBase

if VERSION < v"1.3"
    @warn "To use all the features of Kinetic, please upgrade to Julia 1.3 or newer."
end

export KB

import Base: *
import BSON
import JLD2
import NonlinearSolve
import Roots: Order1, find_zero
import SciMLNLSolve: NLSolveJL
using Base.Threads: @threads
using CSV
using CUDA
using Dates
using Distributed
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
using Printf
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
include("Type/type.jl")
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

const KB = KitBase

end
