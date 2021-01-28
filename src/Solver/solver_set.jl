"""
    struct SolverSet{
        TS<:AbstractSetup,
        TP<:AbstractPhysicalSpace,
        TV<:Union{AbstractVelocitySpace,Nothing},
        TG<:AbstractProperty,
        TI<:AbstractCondition,
        TO<:AbstractString,
    } <: AbstractSolverSet
        # setup
        set::TS
        # physical space
        pSpace::TP
        # velocity space
        vSpace::TV
        # gas property
        gas::TG
        # initial and boundary condition
        ib::TI
        # file system
        outputFolder::TO
    end

Structure of solver setup

"""
struct SolverSet{
    TS<:AbstractSetup,
    TP<:AbstractPhysicalSpace,
    TV<:Union{AbstractVelocitySpace,Nothing},
    TG<:AbstractProperty,
    TI<:AbstractCondition,
    TO<:AbstractString,
} <: AbstractSolverSet
    # setup
    set::TS
    # physical space
    pSpace::TP
    # velocity space
    vSpace::TV
    # gas property
    gas::TG
    # initial and boundary condition
    ib::TI
    # file system
    outputFolder::TO
end

function SolverSet(configfilename::T) where {T<:AbstractString}
    # generate variables from configuration file
    dict = read_dict(configfilename)

    # setup
    set = set_setup(dict)

    # physical space
    pSpace = set_geometry(dict)

    # velocity space
    vSpace = set_velocity(dict)

    # gas property
    gas = set_property(dict)

    # initial & boundary condition
    ib = set_ib(dict, set, vSpace, gas)

    # create working directory
    identifier = string(Dates.now(), "/")
    #outputFolder = string("../out/", replace(identifier, ":"=>"."))
    outputFolder = replace(identifier, ":" => ".")
    mkdir(outputFolder)
    mkdir(string(outputFolder, "data/"))
    cp(configfilename, string(outputFolder, "config.txt"))

    # create new struct
    return SolverSet{
        typeof(set),
        typeof(pSpace),
        typeof(vSpace),
        typeof(gas),
        typeof(ib),
        typeof(outputFolder),
    }(
        set,
        pSpace,
        vSpace,
        gas,
        ib,
        outputFolder,
    )
end

function SolverSet(dict::T) where {T<:AbstractDict}
    # setup
    set = set_setup(;dict...)

    # physical space
    pSpace = set_geometry(dict)

    # velocity space
    vSpace = set_velocity(dict)

    # gas property
    gas = set_property(dict)

    # initial & boundary condition
    ib = set_ib(dict, set, vSpace, gas)

    # create working directory
    identifier = string(Dates.now(), "/")
    outputFolder = replace(identifier, ":" => ".")
    mkdir(outputFolder)
    mkdir(string(outputFolder, "data/"))
    CSV.write(string(outputFolder, "config.csv"), dict)

    # create new struct
    return SolverSet{
        typeof(set),
        typeof(pSpace),
        typeof(vSpace),
        typeof(gas),
        typeof(ib),
        typeof(outputFolder),
    }(
        set,
        pSpace,
        vSpace,
        gas,
        ib,
        outputFolder,
    )
end


"""
    set_setup(dict::T) where {T<:AbstractDict}

Generate AbstractPhysicalSpace

"""
function set_setup(dict::T) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    set = Setup(
        matter,
        case,
        space,
        flux,
        collision,
        nSpecies,
        interpOrder,
        limiter,
        boundary,
        cfl,
        maxTime,
    )

    return set
end

function set_setup(;
    matter,
    case,
    space,
    flux,
    collision,
    nSpecies,
    interpOrder,
    limiter,
    boundary,
    cfl,
    maxTime,
    kwargs...,
)
    set = Setup(
        matter,
        case,
        space,
        flux,
        collision,
        nSpecies,
        interpOrder,
        limiter,
        boundary,
        cfl,
        maxTime,
    )

    return set
end


"""
    set_geometry(dict::T) where {T<:AbstractDict}

Generate AbstractPhysicalSpace

"""
function set_geometry(dict::T) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    Dx = parse(Int, space[1])
    if Dx == 1
        pSpace = PSpace1D(x0, x1, nx, nxg)
    elseif Dx == 2
        pSpace = PSpace2D(x0, x1, nx, y0, y1, ny, nxg, nyg)
    else
        throw("No preset available for 3D simulation, please set it up manually.")
    end

    return pSpace
end

function set_geometry(;
    space,
    x0,
    x1,
    nx,
    nxg,
    y0 = nothing,
    y1 = nothing,
    ny = nothing,
    nyg = nothing,
    kwargs...,
)

    Dx = parse(Int, space[1])
    if Dx == 1
        pSpace = PSpace1D(x0, x1, nx, nxg)
    elseif Dx == 2
        pSpace = PSpace2D(x0, x1, nx, y0, y1, ny, nxg, nyg)
    else
        throw("No preset available for 3D simulation, please set it up manually.")
    end

    return pSpace
end


"""
    set_velocity(dict::T) where {T<:AbstractDict}

Generate AbstractVelocitySpace

"""
function set_velocity(dict::T) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    Dv = parse(Int, space[5])
    if Dv == 0
        vSpace = nothing
    elseif Dv == 1
        if nSpecies == 1
            vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            vSpace = MVSpace1D(umin, umax, ue0, ue1, nu, vMeshType, nug)
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 2
        if nSpecies == 1
            vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            vSpace = MVSpace2D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv,
                vMeshType,
                nug,
                nvg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 3
        if nSpecies == 1
            vSpace = VSpace3D(
                umin,
                umax,
                nu,
                vmin,
                vmax,
                nv,
                wmin,
                wmax,
                nw,
                vMeshType,
                nug,
                nvg,
                nwg,
            )
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            we0 = wmin * sqrt(mi / me)
            we1 = wmax * sqrt(mi / me)
            vSpace = MVSpace3D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv,
                wmin,
                wmax,
                we0,
                we1,
                nw,
                vMeshType,
                nug,
                nvg,
                nwg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    end

    return vSpace
end

function set_velocity(;
    space,
    nSpecies,
    umin,
    umax,
    nu,
    vMeshType,
    nug,
    mi = nothing,
    me = nothing,
    vmin = nothing,
    vmax = nothing,
    nv = nothing,
    nvg = nothing,
    wmin = nothing,
    wmax = nothing,
    nw = nothing,
    nwg = nothing,
    kwargs...,
)
    Dv = parse(Int, space[5])
    if Dv == 0
        vSpace = nothing
    elseif Dv == 1
        if nSpecies == 1
            vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            vSpace = MVSpace1D(umin, umax, ue0, ue1, nu, vMeshType, nug)
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 2
        if nSpecies == 1
            vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            vSpace = MVSpace2D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv,
                vMeshType,
                nug,
                nvg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 3
        if nSpecies == 1
            vSpace = VSpace3D(
                umin,
                umax,
                nu,
                vmin,
                vmax,
                nv,
                wmin,
                wmax,
                nw,
                vMeshType,
                nug,
                nvg,
                nwg,
            )
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            we0 = wmin * sqrt(mi / me)
            we1 = wmax * sqrt(mi / me)
            vSpace = MVSpace3D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv,
                wmin,
                wmax,
                we0,
                we1,
                nw,
                vMeshType,
                nug,
                nvg,
                nwg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    end

    return vSpace
end


"""
    set_property(dict::T) where {T<:AbstractDict}

Generate AbstractProperty

"""
function set_property(dict::T) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    Dx = parse(Int, space[1])
    γD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
        if x == 0
            return Dx
        elseif x == 1 # 1f
            if y >= 3 # 3v
                return 3
            else # 1v
                return Dx
            end
        elseif x == 2 # 2f
            return Dx
        elseif x >= 3 # 3f / 4f
            return 3
        else
            return nothing
        end
    end
    γ = heat_capacity_ratio(inK, γD)

    if matter == "gas"

        if nSpecies == 1
            μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)
            gas = Gas(knudsen, mach, prandtl, inK, γ, omega, alphaRef, omegaRef, μᵣ)
        elseif nSpecies == 2
            kne = knudsen * (me / mi)
            gas = Mixture([knudsen, kne], mach, prandtl, inK, γ, mi, ni, me, ne)
        else
            throw("The gas property only supports up to two species.")
        end
    
    elseif matter == "plasma"

        kne = knudsen * (me / mi)
        if Dx == 1
            gas = Plasma1D(
                [knudsen, kne],
                mach,
                prandtl,
                inK,
                γ,
                mi,
                ni,
                me,
                ne,
                lD,
                rL,
                sol,
                echi,
                bnu,
            )
        elseif Dx == 2
            gas = Plasma2D(
                [knudsen, kne],
                mach,
                prandtl,
                inK,
                γ,
                mi,
                ni,
                me,
                ne,
                lD,
                rL,
                sol,
                echi,
                bnu,
            )
        else
            throw("The plasma property only supports up to 2D case.")
        end
    
    end

    return gas
end


"""
    set_ib(dict::T, set, vSpace, gas) where {T<:AbstractDict}

    set_ib(
        set::T,
        vSpace,
        gas,
        uLid = 0.15,
        vLid = 0.0,
        tLid = 1.0,
    ) where {T<:AbstractSetup}

Generate AbstractIB

"""
function set_ib(dict::T, set, vSpace, gas) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    if @isdefined uLid
        ib = set_ib(set, vSpace, gas, uLid, vLid, tLid)
    else
        ib = set_ib(set, vSpace, gas)
    end

    return ib
end

function set_ib(
    set::T,
    vSpace,
    gas,
    uLid = 0.15,
    vLid = 0.0,
    tLid = 1.0,
) where {T<:AbstractSetup}

    if set.case == "shock"

        if set.nSpecies == 1
            if set.space[3:end] == "1f1v"
                wL, primL, fL, bcL, wR, primR, fR, bcR = ib_rh(gas.Ma, gas.γ, vSpace.u)
                ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
            elseif set.space[3:end] == "1f3v"
                wL, primL, fL, bcL, wR, primR, fR, bcR =
                    ib_rh(gas.Ma, gas.γ, vSpace.u, vSpace.v, vSpace.w)
                ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
            elseif set.space[3:end] == "2f1v"
                wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                    ib_rh(gas.Ma, gas.γ, gas.K, vSpace.u)
                ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
            end
        elseif set.nSpecies == 2
            if set.space[3:end] == "2f1v"
                wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                    ib_rh(gas.Ma, gas.γ, gas.K, gas.mi, gas.me, gas.ni, gas.ne, vSpace.u)
                ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
            end
        end

    elseif set.case == "sod"

        if set.nSpecies == 1
            if set.space[3:end] == "0f0v"
                wL, primL, bcL, wR, primR, bcR = ib_sod(gas.γ)
                ib = IB(wL, primL, bcL, wR, primR, bcR)
            elseif set.space[3:end] == "1f1v"
                wL, primL, fL, bcL, wR, primR, fR, bcR = ib_sod(gas.γ, vSpace.u)
                ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
            elseif set.space[3:end] == "1f3v"
                wL, primL, fL, bcL, wR, primR, fR, bcR =
                    ib_sod(gas.γ, vSpace.u, vSpace.v, vSpace.w)
                ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
            elseif set.space[3:end] == "2f1v"
                wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR = ib_sod(gas.γ, gas.K, vSpace.u)
                ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
            end
        elseif set.nSpecies == 2
            if set.space[3:end] == "0f0v"
                wL, primL, bcL, wR, primR, bcR = ib_sod(gas.γ, gas.mi, gas.me)
                ib = IB(wL, primL, bcL, wR, primR, bcR)
            end
        end

    elseif set.case == "brio-wu"

        if set.space[3:end] == "4f1v"
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
            lorenzR = ib_briowu(gas.γ, gas.mi, gas.me, vSpace.u)

            ib = IB4F(
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
            )
        elseif set.space[3:end] == "3f2v"
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
            lorenzR = ib_briowu(gas.γ, gas.mi, gas.me, vSpace.u, vSpace.v)

            ib = IB3F(
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
            )
        end

    elseif set.case == "cavity"

        if set.space[3:end] == "1f2v"
            wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD =
                ib_cavity(gas.γ, uLid, vLid, tLid, vSpace.u, vSpace.v)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD)
        elseif set.space[3:end] == "2f2v"
            wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD =
                ib_cavity(gas.γ, gas.K, uLid, vLid, tLid, vSpace.u, vSpace.v)
            ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD)
        end

    else

        throw("No default ib available, please set it up manually.")

    end

    return ib

end
