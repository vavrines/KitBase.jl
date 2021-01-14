"""
KitBase.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2021 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitBase

if VERSION < v"1.3"
    @warn "Kinetic.jl matches perfectly with Julia 1.3 or newer versions."
end

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
using Optim
using MultivariatePolynomials
using TypedPolynomials
using PyCall
using CUDA

include("Data/data.jl")
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
    threads = Threads.nthreads()
    if threads == 1
        @info "Kinetic will run serially"
    elseif threads > 1
        @info "Kinetic will run with $threads threads"
        # https://github.com/CliMA/Oceananigans.jl/issues/1113
        FFTW.set_num_threads(4 * threads)
    end

    if has_cuda()
        @info "Kinetic will enable CUDA devices"
        for (i, dev) in enumerate(CUDA.devices())
            @info "$i: $(CUDA.name(dev))"
        end

        @info "Scalar operation is disabled in CUDA"
        CUDA.allowscalar(false)
    end
end

end
