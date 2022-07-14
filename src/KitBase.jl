"""
KitBase.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2022 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitBase

if VERSION < v"1.3"
    @warn "To use all the features of Kinetic, please upgrade to Julia 1.3 or newer."
end

export KB

import Base: *, @kwdef
import BSON
import JLD2
using Base.Threads: @threads
using CSV
using CUDA
using Dates
using DiffEqOperators
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
using Printf
using PyCall
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
const itp = PyNULL()

function __init__()
    np = nworkers()
    nt = Threads.nthreads()
    if np == 1 && nt == 1 && !has_cuda()
        @info "Kinetic will run in serial"
    else
        if has_cuda()
            @info "Kinetic will run in parallel with $np workers, $nt threads and CUDA"
            for (i, dev) in enumerate(CUDA.devices())
                @info "$i: $(CUDA.name(dev))"
            end
            #@info "Scalar operation is disabled in CUDA"
            CUDA.allowscalar(false)
        else
            @info "Kinetic will run in parallel with $np workers and $nt threads"
        end
    end

    copy!(itp, pyimport("scipy.interpolate"))
end

end
