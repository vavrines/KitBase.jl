# ============================================================
# KitBase.jl: The lightweight prototype with physical 
#             formulations of Kinetic.jl Ecosystem
# Copyright (c) 2020 Tianbai Xiao <tianbaixiao@gmail.com>
# ============================================================

module KitBase

using Dates
using OffsetArrays
using LinearAlgebra
using FastGaussQuadrature
using SpecialFunctions
using FFTW
using Plots
using FileIO
using JLD2
using ProgressMeter
using Distributions
using PyCall

include("Type/type.jl")
include("IO/io.jl")
include("Math/math.jl")
include("Geometry/geometry.jl")
include("Theory/theory.jl")
include("Phase/phase.jl")
include("Reconstruction/reconstruction.jl")
include("Flux/flux.jl")
include("Config/config.jl")
include("Solver/solver.jl")

end
