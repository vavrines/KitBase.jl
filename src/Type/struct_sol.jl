# ============================================================
# Structs of Solution
# Solver stores data in structs of arrays (SoA)
# ============================================================

struct Solution{T1,T2,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
end

struct Solution1F{T1,T2,T3,T4,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
    f::T3
    ∇f::T4
end

struct Solution2F{T1,T2,T3,T4,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
    h::T3
    b::T3
    ∇h::T4
    ∇b::T4
end


function Solution1D(w::AA, prim::AA)
    sw = zero(w)
    return Solution{typeof(w),typeof(sw),1}(w, prim, sw)
end

function Solution1D(w::AA, prim::AA, f::AA)
    sw = zero(w)
    sf = zero(f)
    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),1}(w, prim, sw, f, sf)
end

function Solution1D(w::AA, prim::AA, h::AA, b::AA)
    sw = zero(w)
    sh = zero(h)
    sb = zero(b)
    return Solution2F{typeof(w),typeof(sw),typeof(h),typeof(sh),1}(w, prim, sw, h, b, sh, sb)
end


function Solution2D(w::AA, prim::AA)
    sw = zeros(eltype(w), axes(w)[1], 2, axes(w)[2:end]...)
    return Solution{typeof(w),typeof(sw),2}(w, prim, sw)
end

function Solution2D(w::AA, prim::AA, f::AA)
    sw = zeros(eltype(w), axes(w)[1], 2, axes(w)[2:end]...)
    sf = begin
        if ndims(f) - ndims(w) > 1
            zeros(eltype(f), axes(f)[1:2]..., 2, axes(f)[3:end]...)
        else
            zeros(eltype(f), axes(f)[1], 2, axes(f)[2:end]...)
        end
    end
    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),2}(w, prim, sw, f, sf)
end

function Solution2D(w::AA, prim::AA, h::AA, b::AA)
    sw = zeros(eltype(w), axes(w)[1], 2, axes(w)[2:end]...)
    sh = begin
        if ndims(h) - ndims(w) > 1
            zeros(eltype(h), axes(h)[1:2]..., 2, axes(h)[3:end]...)
        else
            zeros(eltype(h), axes(h)[1], 2, axes(h)[2:end]...)
        end
    end
    sb = zero(sh)
    return Solution2F{typeof(w),typeof(sw),typeof(h),typeof(sh),2}(w, prim, sw, h, b, sh, sb)
end
