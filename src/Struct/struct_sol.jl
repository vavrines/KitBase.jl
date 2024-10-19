# ============================================================
# Structs of Solution
# Solver stores data in structs of arrays (SoA)
# ============================================================

"""
$(TYPEDEF)

Solution structure with no distribution function
"""
struct Solution{T1,T2,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
end

"""
$(TYPEDEF)

Solution structure with 1 distribution function
"""
struct Solution1F{T1,T2,T3,T4,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
    f::T3
    ∇f::T4
end

"""
$(TYPEDEF)

Solution structure with 2 distribution functions
"""
struct Solution2F{T1,T2,T3,T4,ND} <: AbstractSolution
    w::T1
    prim::T1
    ∇w::T2
    h::T3
    b::T3
    ∇h::T4
    ∇b::T4
end

"""
$(SIGNATURES)

Generate 1D solution structure
"""
function Solution1D(w::AA, prim::AA)
    sw = @. zero(w)
    return Solution{typeof(w),typeof(sw),1}(w, prim, sw)
end

"""
$(SIGNATURES)
"""
function Solution1D(w::AA, prim::AA, f::AA)
    sw = @. zero(w)
    sf = @. zero(f)

    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),1}(w, prim, sw, f, sf)
end

"""
$(SIGNATURES)
"""
function Solution1D(w::AA, prim::AA, h::AA, b::AA)
    sw = @. zero(w)
    sh = @. zero(h)
    sb = @. zero(b)

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

"""
$(SIGNATURES)

Generate 2D solution structure
"""
function Solution2D(w::AA, prim::AA)
    if w[1] isa Number
        sw = slope_array(w)
    else
        if ndims(w) == 1
            sw = [slope_array(w[i]) for i in axes(w, 1)]
        else
            sw = [slope_array(w[i, j]) for i in axes(w, 1), j in axes(w, 2)]
        end
    end

    return Solution{typeof(w),typeof(sw),2}(w, prim, sw)
end

"""
$(SIGNATURES)
"""
function Solution2D(w::AA, prim::AA, f::AA)
    if w[1] isa Number
        sw = slope_array(w)
        sf = slope_array(f)
    else
        if ndims(w) == 1
            sw = [slope_array(w[i]) for i in axes(w, 1)]
            sf = [slope_array(f[i]) for i in axes(w, 1)]
        else
            sw = [slope_array(w[i, j]) for i in axes(w, 1), j in axes(w, 2)]
            sf = [slope_array(f[i, j]) for i in axes(w, 1), j in axes(w, 2)]
        end
    end

    return Solution1F{typeof(w),typeof(sw),typeof(f),typeof(sf),2}(w, prim, sw, f, sf)
end

"""
$(SIGNATURES)
"""
function Solution2D(w::AA, prim::AA, h::AA, b::AA)
    if w[1] isa Number
        sw = slope_array(w)
        sh = slope_array(h)
        sb = slope_array(b)
    else
        if ndims(w) == 1
            sw = [slope_array(w[i]) for i in axes(w, 1)]
            sh = [slope_array(h[i]) for i in axes(w, 1)]
            sb = [slope_array(b[i]) for i in axes(w, 1)]
        else
            sw = [slope_array(w[i, j]) for i in axes(w, 1), j in axes(w, 2)]
            sh = [slope_array(h[i, j]) for i in axes(w, 1), j in axes(w, 2)]
            sb = [slope_array(b[i, j]) for i in axes(w, 1), j in axes(w, 2)]
        end
    end

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
