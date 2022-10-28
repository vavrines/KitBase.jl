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

function SolverSet(file)
    # generate configuration
    config = begin
        if file isa AbstractString
            read_cfg(file)
        else
            file
        end
    end

    # setup
    set = set_setup(config)

    # physical space
    ps = set_geometry(config)

    # velocity space
    vs = set_velocity(config)

    # gas property
    gas = set_property(config)

    # initial & boundary condition
    ib = set_ib(config, set, ps, vs, gas)

    # create working directory
    identifier = string(Dates.now(), "/")
    #outputFolder = string("../out/", replace(identifier, ":"=>"."))
    outputFolder = replace(identifier, ":" => ".")
    mkdir(outputFolder)
    mkdir(string(outputFolder, "data/"))
    if file isa AbstractString
        cp(configfilename, string(outputFolder, "config.txt"))
    else
        if file isa NamedTuple
            dict = ntuple_dict(file)
        else
            dict = file
        end
        CSV.write(string(outputFolder, "config.csv"), dict)
    end

    # create new struct
    return SolverSet(set, ps, vs, gas, ib, outputFolder)
end


"""
$(SIGNATURES)

Generate AbstractPhysicalSpace
"""
set_setup(dict::Union{AbstractDict,NamedTuple}) = set_setup(; dict...)

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
set_geometry(dict::Union{AbstractDict,NamedTuple}) = set_geometry(; dict...)

"""
$(SIGNATURES)
"""
function set_geometry(;
    space,
    x0 = nothing,
    x1 = nothing,
    nx = nothing,
    nxg = nothing,
    y0 = nothing,
    y1 = nothing,
    ny = nothing,
    nyg = nothing,
    mesh = nothing,
    kwargs...,
)
    try
        return UnstructPSpace(mesh)
    catch
        Dx = parse(Int, space[1])
        if Dx == 1
            return PSpace1D(x0, x1, nx, nxg)
        elseif Dx == 2
            return PSpace2D(x0, x1, nx, y0, y1, ny, nxg, nyg)
        else
            throw("No preset available for 3D simulation, please set it up manually.")
        end
    end
end


"""
$(SIGNATURES)

Generate AbstractVelocitySpace
"""
set_velocity(dict::Union{AbstractDict,NamedTuple}) = set_velocity(; dict...)

"""
$(SIGNATURES)
"""
function set_velocity(;
    space,
    nSpecies,
    umin = nothing,
    umax = nothing,
    nu = nothing,
    vMeshType = nothing,
    nug = nothing,
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
            vSpace = VSpace1D(umin, umax, nu; type = vMeshType, ng = nug)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            vSpace = MVSpace1D(umin, umax, ue0, ue1, nu; type = vMeshType, ng = nug)
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 2
        if nSpecies == 1
            vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv; type = vMeshType, ngu = nug, ngv = nvg)
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
                nv;
                type = vMeshType,
                ngu = nug,
                ngv = nvg,
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
                nw;
                type = vMeshType,
                ngu = nug,
                ngv = nvg,
                ngw = nwg,
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
                nw;
                type = vMeshType,
                ngu = nug,
                ngv = nvg,
                ngw = nwg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    end

    return vSpace
end


"""
$(SIGNATURES)

Generate property of matter
"""
function set_property(dict::Union{AbstractDict,NamedTuple})
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

        if dict[:nSpecies] == 1
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
        elseif dict[:nSpecies] == 2
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

    elseif dict[:matter] == "plasma"

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

Generate initial & boundary conditions
"""
function set_ib(dict::Union{AbstractDict,NamedTuple}, set, ps, vs, gas)
    if haskey(dict, :uLid)
        ib = set_ib(set, ps, vs, gas, dict[:uLid], dict[:vLid], dict[:tLid])
    else
        ib = set_ib(set, ps, vs, gas)
    end

    return ib
end

"""
$(SIGNATURES)
"""
function set_ib(set, pSpace, vSpace, gas, Um = 0.15, Vm = 0.0, Tm = 1.0)
    ib = begin
        if parse(Int, set.space[3]) == 0
            fw, bc, p = config_ib(set, pSpace, vSpace, gas)
            IB{typeof(bc)}(fw, bc, p)
        elseif parse(Int, set.space[3]) in [3, 4] && gas isa AbstractPlasma
            fw, ff, fE, fB, fL, bc, p = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname)){typeof(bc)}(fw, ff, fE, fB, fL, bc, p)
        elseif set.case == "cavity"
            fw, ff, bc, p = config_ib(set, pSpace, vSpace, gas, Um, Vm, Tm)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname)){typeof(bc)}(fw, ff, bc, p)
        else
            fw, ff, bc, p = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname)){typeof(bc)}(fw, ff, bc, p)
        end
    end

    return ib
end
