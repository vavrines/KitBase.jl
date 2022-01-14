"""
$(TYPEDEF)

Structure of solver setup

# Fields

$(FIELDS)
"""
struct SolverSet{
    TS<:AbstractSetup,
    TP<:AbstractPhysicalSpace,
    TV<:Union{AbstractVelocitySpace,Nothing},
    TG<:AbstractProperty,
    TI<:Union{AbstractCondition,Nothing},
    TO<:AbstractString,
} <: AbstractSolverSet
    # setup
    set::TS
    # physical space
    pSpace::TP
    ps::TP
    # velocity space
    vSpace::TV
    vs::TV
    # gas property
    gas::TG
    # initial and boundary condition
    ib::TI
    # file system
    outputFolder::TO
end

function SolverSet(
    SET::AbstractSetup,
    PS::AbstractPhysicalSpace,
    VS::Union{AbstractVelocitySpace,Nothing},
    GAS::AbstractProperty,
    IB::Union{AbstractCondition,Nothing},
    DIR::AbstractString = @__DIR__,
)
    return SolverSet{typeof(SET),typeof(PS),typeof(VS),typeof(GAS),typeof(IB),typeof(DIR)}(
        SET,
        PS,
        PS,
        VS,
        VS,
        GAS,
        IB,
        DIR,
    )
end

function SolverSet(configfilename::T) where {T<:AbstractString}
    # generate variables from configuration file
    dict = read_cfg(configfilename)

    # setup
    set = set_setup(dict)

    # physical space
    ps = set_geometry(dict)

    # velocity space
    vs = set_velocity(dict)

    # gas property
    gas = set_property(dict)

    # initial & boundary condition
    ib = set_ib(dict, set, ps, vs, gas)

    # create working directory
    identifier = string(Dates.now(), "/")
    #outputFolder = string("../out/", replace(identifier, ":"=>"."))
    outputFolder = replace(identifier, ":" => ".")
    mkdir(outputFolder)
    mkdir(string(outputFolder, "data/"))
    cp(configfilename, string(outputFolder, "config.txt"))

    # create new struct
    return SolverSet(set, ps, vs, gas, ib, outputFolder)
end

function SolverSet(dict::T) where {T<:AbstractDict}
    # setup
    set = set_setup(; dict...)

    # physical space
    ps = set_geometry(; dict...)

    # velocity space
    vs = set_velocity(; dict...)

    # gas property
    gas = set_property(dict)

    # initial & boundary condition
    ib = set_ib(dict, set, ps, vs, gas)

    # create working directory
    identifier = string(Dates.now(), "/")
    outputFolder = replace(identifier, ":" => ".")
    mkdir(outputFolder)
    mkdir(string(outputFolder, "data/"))
    CSV.write(string(outputFolder, "config.csv"), dict)

    # create new struct
    return SolverSet(set, ps, vs, gas, ib, outputFolder)
end


"""
$(SIGNATURES)

Generate AbstractPhysicalSpace
"""
function set_setup(dict::T) where {T<:AbstractDict}
    set = Setup(
        dict[:matter],
        dict[:case],
        dict[:space],
        dict[:flux],
        dict[:collision],
        dict[:nSpecies],
        dict[:interpOrder],
        dict[:limiter],
        dict[:boundary],
        dict[:cfl],
        dict[:maxTime],
    )

    return set
end

"""
$(SIGNATURES)
"""
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
$(SIGNATURES)

Generate AbstractPhysicalSpace
"""
function set_geometry(dict::T) where {T<:AbstractDict}
    try
        return UnstructPSpace(dict[:mesh])
    catch
        if parse(Int, dict[:space][1]) == 1
            return PSpace1D(dict[:x0], dict[:x1], dict[:nx], dict[:nxg])
        elseif parse(Int, dict[:space][1]) == 2
            return PSpace2D(
                dict[:x0],
                dict[:x1],
                dict[:nx],
                dict[:y0],
                dict[:y1],
                dict[:ny],
                dict[:nxg],
                dict[:nyg],
            )
        else
            throw("No preset available for 3D simulation, please set it up manually.")
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
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
$(SIGNATURES)

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

"""
$(SIGNATURES)
"""
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
$(SIGNATURES)

Generate AbstractProperty
"""
function set_property(dict::T) where {T<:AbstractDict}
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end

    Dx = parse(Int, dict[:space][1])
    γD = map(parse(Int, dict[:space][3]), parse(Int, dict[:space][5])) do x, y # (x)f(y)v
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
    γ = heat_capacity_ratio(dict[:inK], γD)

    if dict[:matter] == "radiation"

        gas = Radiation(dict[:knudsen], dict[:sigmaS], dict[:sigmaA])

    elseif dict[:matter] == "gas"

        if nSpecies == 1
            μᵣ = ref_vhs_vis(dict[:knudsen], dict[:alphaRef], dict[:omegaRef])

            if dict[:collision] == "fsm"
                nm = begin
                    if haskey(dict, :nm)
                        dict[:nm]
                    else
                        return 5
                    end
                end
                vs = set_velocity(; dict...)
                fsm = fsm_kernel(vs, μᵣ, nm, dict[:alphaRef])
            else
                fsm = nothing
            end

            gas = Gas(
                Kn = dict[:knudsen],
                Ma = dict[:mach],
                Pr = dict[:prandtl],
                K = dict[:inK],
                γ = γ,
                ω = dict[:omega],
                αᵣ = dict[:alphaRef],
                ωᵣ = dict[:omegaRef],
                μᵣ = μᵣ,
                fsm = fsm,
            )
        elseif nSpecies == 2
            kne = dict[:knudsen] * (dict[:me] / dict[:mi])
            gas = Mixture(
                [dict[:knudsen], kne],
                dict[:mach],
                dict[:prandtl],
                dict[:inK],
                γ,
                dict[:mi],
                dict[:ni],
                dict[:me],
                dict[:ne],
            )
        else
            throw("The gas property only supports up to two species.")
        end

    elseif matter == "plasma"

        kne = dict[:knudsen] * (dict[:me] / dict[:mi])
        if Dx == 1
            gas = Plasma1D(
                [dict[:knudsen], kne],
                dict[:mach],
                dict[:prandtl],
                dict[:inK],
                γ,
                dict[:mi],
                dict[:ni],
                dict[:me],
                dict[:ne],
                dict[:lD],
                dict[:rL],
                dict[:sol],
                dict[:echi],
                dict[:bnu],
            )
        elseif Dx == 2
            gas = Plasma2D(
                [dict[:knudsen], kne],
                dict[:mach],
                dict[:prandtl],
                dict[:inK],
                γ,
                dict[:mi],
                dict[:ni],
                dict[:me],
                dict[:ne],
                dict[:lD],
                dict[:rL],
                dict[:sol],
                dict[:echi],
                dict[:bnu],
            )
        else
            throw("The plasma property only supports up to 2D case.")
        end

    end

    return gas
end


"""
$(SIGNATURES)

Generate AbstractIB
"""
function set_ib(dict::T, set, ps, vs, gas) where {T<:AbstractDict}
    if haskey(dict, :uLid)
        ib = set_ib(set, ps, vs, gas, uLid, vLid, tLid)
    else
        ib = set_ib(set, ps, vs, gas)
    end

    return ib
end

"""
$(SIGNATURES)
"""
function set_ib(set, pSpace, vSpace, gas, Um = 0.15, Vm = 0.0, Tm = 1.0)
    fname = begin
        if set.case == "shock"
            "ib_rh"
        elseif set.case == "brio-wu"
            "ib_briowu"
        else
            "ib_" * string(set.case)
        end
    end

    ib = begin
        if parse(Int, set.space[3]) == 0
            fw, bc = config_ib(set, pSpace, vSpace, gas)
            IB(fw, bc)
        elseif parse(Int, set.space[3]) in [3, 4] && gas isa AbstractPlasma
            fw, ff, fE, fB, fL, bc = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, fE, fB, fL, bc)
        elseif set.case == "cavity"
            fw, ff, bc = config_ib(set, pSpace, vSpace, gas, Um, Vm, Tm)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, bc)
        else
            fw, ff, bc = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, bc)
        end
    end

    return ib
end
