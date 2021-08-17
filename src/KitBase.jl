"""
KitBase.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2021 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitBase

if VERSION < v"1.3"
    @warn "Kinetic works better with Julia 1.3 or newer versions."
end

using Base.Threads: @threads
using Reexport
@reexport using FiniteMesh
using CSV
using CUDA
using Dates
using Distributed
using Distributions
using FastGaussQuadrature
using FFTW
using FileIO
using JLD2
using LinearAlgebra
using MultivariatePolynomials
using OffsetArrays
using Optim
using Plots
using ProgressMeter
using PyCall
using SpecialFunctions
using StaticArrays
using StructArrays
using TypedPolynomials
using WriteVTK

const itp = PyNULL()

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
include("Config/config.jl")
include("Solver/solver.jl")

function __init__()
    np = nworkers()
    nt = Threads.nthreads()
    if nt > 1 || np > 1
        @info "Kinetic will run with $np processors and $nt threads"
    else
        @info "Kinetic will run serially"
    end

    if has_cuda()
        @info "Kinetic will run with CUDA"
        for (i, dev) in enumerate(CUDA.devices())
            @info "$i: $(CUDA.name(dev))"
        end
        @info "Scalar operation is disabled in CUDA"
        CUDA.allowscalar(false)
    end

    copy!(itp, pyimport("scipy.interpolate"))
end

end
