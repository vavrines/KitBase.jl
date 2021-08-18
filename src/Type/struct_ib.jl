# ============================================================
# Structs of Initial & Boundary Conditions
# ============================================================

"""
    mutable struct IB{A,B,C} <: AbstractCondition
        wL::A
        primL::B
        bcL::B

        wR::A
        primR::B
        bcR::B

        bcU::B
        bcD::B

        vL::C
        vR::C
    end

Initial & boundary condition with no distribution function

"""
mutable struct IB{A,B,C} <: AbstractCondition
    wL::A
    primL::B
    bcL::B
    wR::A
    primR::B
    bcR::B
    bcU::B
    bcD::B
    vL::C
    vR::C
end

# works for both 1V/3V and single-/multi-component gases
function IB(wL, primL, bcL, wR, primR, bcR, bcU = deepcopy(bcR), bcD = deepcopy(bcR))

    if ndims(primL) == 0
        vL = primL
        vR = primR
    elseif ndims(primL) == 1
        vL = zeros(eltype(primL), 3)
        vR = zeros(eltype(primR), 3)

        if size(primL, 1) == 3
            vL[1] = primL[2]
            vR[1] = primR[2]
        elseif size(primL, 1) == 4
            vL[1:2] .= primL[2:3]
            vR[1:2] .= primR[2:3]
        elseif size(primL, 1) == 5
            vL .= primL[2:4]
            vR .= primR[2:4]
        end
    elseif ndims(primL) == 2
        vL = zeros(eltype(primL), 3, size(primL, 2))
        vR = zeros(eltype(primR), 3, size(primR, 2))

        if size(primL, 1) == 3
            vL[1, :] .= primL[2, :]
            vR[1, :] .= primR[2, :]
        elseif size(primL, 1) == 4
            vL[1:2, :] .= primL[2:3, :]
            vR[1:2, :] .= primR[2:3, :]
        elseif size(primL, 1) == 5
            vL .= primL[2:4, :]
            vR .= primR[2:4, :]
        end
    end

    return IB(wL, primL, bcL, wR, primR, bcR, bcU, bcD, vL, vR)
end


"""
    mutable struct IB1F{A,B} <: AbstractCondition
        wL::A
        primL::A
        fL::B
        bcL::A

        wR::A
        primR::A
        fR::B
        bcR::A

        bcU::A
        bcD::A
    end

Initial & boundary condition with 1 distribution function

"""
mutable struct IB1F{A<:Function,B} <: AbstractCondition
    fw::A
    ff::A
    bcL::B
    bcR::B
    bcU::B
    bcD::B
end

# works for both 1V/3V and single-/multi-component gases
function IB1F(
    fw::Function,
    ff::Function,
    bcL,
    bcR,
)
    return IB1F{Function,typeof(bcL)}(fw, ff, bcL, bcR, bcR, bcR)
end


"""
    mutable struct IB2F{A,B} <: AbstractCondition
        wL::A
        primL::A
        hL::B
        bL::B
        bcL::A

        wR::A
        primR::A
        hR::B
        bR::B
        bcR::A

        bcU::A
        bcD::A
    end

Initial & boundary condition with 2 distribution functions

"""
mutable struct IB2F{A,B} <: AbstractCondition
    wL::A
    primL::A
    hL::B
    bL::B
    bcL::A
    wR::A
    primR::A
    hR::B
    bR::B
    bcR::A
    bcU::A
    bcD::A

    function IB2F(
        wL,
        primL,
        hL,
        bL,
        bcL,
        wR,
        primR,
        hR,
        bR,
        bcR,
        bcU = deepcopy(bcR),
        bcD = deepcopy(bcR),
    )
        new{typeof(wL),typeof(hL)}(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD)
    end
end


"""
    mutable struct IB3F{A,B,C,D} <: AbstractCondition
        wL::A
        primL::A
        h0L::B
        h1L::B
        h2L::B
        bcL::A
        EL::C
        BL::C
        lorenzL::D

        wR::A
        primR::A
        h0R::B
        h1R::B
        h2R::B
        bcR::A
        ER::C
        BR::C
        lorenzR::D

        bcU::A
        bcD::A
    end

Initial & boundary condition with 3 distribution functions

"""
mutable struct IB3F{A,B,C,D} <: AbstractCondition
    wL::A
    primL::A
    h0L::B
    h1L::B
    h2L::B
    bcL::A
    EL::C
    BL::C
    lorenzL::D

    wR::A
    primR::A
    h0R::B
    h1R::B
    h2R::B
    bcR::A
    ER::C
    BR::C
    lorenzR::D

    bcU::A
    bcD::A

    function IB3F(
        wL::AbstractArray,
        primL::AbstractArray,
        h0L::AbstractArray,
        h1L::AbstractArray,
        h2L::AbstractArray,
        bcL::AbstractArray,
        EL,
        BL,
        lorenzL,
        wR::AbstractArray,
        primR::AbstractArray,
        h0R::AbstractArray,
        h1R::AbstractArray,
        h2R::AbstractArray,
        bcR::AbstractArray,
        ER,
        BR,
        lorenzR,
        bcU = deepcopy(bcR),
        bcD = deepcopy(bcR),
    )
        new{typeof(wL),typeof(h0L),typeof(EL),typeof(lorenzL)}(
            wL,
            primL,
            h0L,
            h1L,
            h2L,
            bcL,
            EL,
            BL,
            lorenzL,
            wR,
            primR,
            h0R,
            h1R,
            h2R,
            bcR,
            ER,
            BR,
            lorenzR,
            bcU,
            bcD,
        )
    end

end


"""
    mutable struct IB4F{A,B,C,D} <: AbstractCondition
        wL::A
        primL::A
        h0L::B
        h1L::B
        h2L::B
        h3L::B
        bcL::A
        EL::C
        BL::C
        lorenzL::D

        wR::A
        primR::A
        h0R::B
        h1R::B
        h2R::B
        h3R::B
        bcR::A
        ER::C
        BR::C
        lorenzR::D

        bcU::A
        bcD::A
    end

Initial & boundary condition with 4 distribution functions

"""
mutable struct IB4F{A,B,C,D} <: AbstractCondition

    wL::A
    primL::A
    h0L::B
    h1L::B
    h2L::B
    h3L::B
    bcL::A
    EL::C
    BL::C
    lorenzL::D

    wR::A
    primR::A
    h0R::B
    h1R::B
    h2R::B
    h3R::B
    bcR::A
    ER::C
    BR::C
    lorenzR::D

    bcU::A
    bcD::A

    function IB4F(
        wL::AbstractArray,
        primL::AbstractArray,
        h0L::AbstractArray,
        h1L::AbstractArray,
        h2L::AbstractArray,
        h3L::AbstractArray,
        bcL::AbstractArray,
        EL::AbstractArray,
        BL::AbstractArray,
        lorenzL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        h0R::AbstractArray,
        h1R::AbstractArray,
        h2R::AbstractArray,
        h3R::AbstractArray,
        bcR::AbstractArray,
        ER::AbstractArray,
        BR::AbstractArray,
        lorenzR::AbstractArray,
        bcU = deepcopy(bcR)::AbstractArray,
        bcD = deepcopy(bcR)::AbstractArray,
    )

        new{typeof(wL),typeof(h0L),typeof(EL),typeof(lorenzL)}(
            wL,
            primL,
            h0L,
            h1L,
            h2L,
            h3L,
            bcL,
            EL,
            BL,
            lorenzL,
            wR,
            primR,
            h0R,
            h1R,
            h2R,
            h3R,
            bcR,
            ER,
            BR,
            lorenzR,
            bcU,
            bcD,
        )

    end

end
