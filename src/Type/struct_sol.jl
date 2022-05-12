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
    sw = deepcopy(w)
    if w[1] isa Number
        sw .= 0.0
    else
        for e in sw
            e .= 0.0
        end
    end

    return Solution{typeof(w),typeof(sw),1}(w, prim, sw)
end

function Solution1D(w::AA, prim::AA, f::AA)
    sw = deepcopy(w)
    sf = deepcopy(f)
    if w[1] isa Number
        sw .= 0.0
        sf .= 0.0
    else
        for e in sw
            e .= 0.0
        end
        for e in sf
            e .= 0.0
        end
    end

    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),1}(w, prim, sw, f, sf)
end

function Solution1D(w::AA, prim::AA, h::AA, b::AA)
    sw = deepcopy(w)
    sh = deepcopy(h)
    sb = deepcopy(b)
    if w[1] isa Number
        sw .= 0.0
        sh .= 0.0
        sb .= 0.0
    else
        for e in sw
            e .= 0.0
        end
        for e in sh
            e .= 0.0
        end
        for e in sb
            e .= 0.0
        end
    end

    return Solution2F{typeof(w),typeof(sw),typeof(h),typeof(sh),1}(
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )
end


function Solution2D(w::AA, prim::AA)
    ET = ifelse(w[1] isa Number, eltype(w), eltype(w[1]))
    sw = zeros(ET, axes(w)[1], 2, axes(w)[2:end]...)
    return Solution{typeof(w),typeof(sw),2}(w, prim, sw)
end

function Solution2D(w::AA, prim::AA, f::AA)
    ET = ifelse(w[1] isa Number, eltype(w), eltype(w[1]))
    sw = zeros(ET, axes(w)[1], 2, axes(w)[2:end]...)
    sf = begin
        if ndims(f) - ndims(w) > 1
            zeros(ET, axes(f)[1:2]..., 2, axes(f)[3:end]...)
        else
            zeros(ET, axes(f)[1], 2, axes(f)[2:end]...)
        end
    end
    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),2}(w, prim, sw, f, sf)
end

function Solution2D(w::AA, prim::AA, h::AA, b::AA)
    ET = ifelse(w[1] isa Number, eltype(w), eltype(w[1]))
    sw = zeros(ET, axes(w)[1], 2, axes(w)[2:end]...)
    sh = begin
        if ndims(h) - ndims(w) > 1
            zeros(ET, axes(h)[1:2]..., 2, axes(h)[3:end]...)
        else
            zeros(ET, axes(h)[1], 2, axes(h)[2:end]...)
        end
    end
    sb = zero(sh)
    return Solution2F{typeof(w),typeof(sw),typeof(h),typeof(sh),2}(
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )
end
