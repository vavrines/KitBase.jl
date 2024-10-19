"""
$(TYPEDEF)

Structure of solver setup

## Fields

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
    IB::Union{AbstractCondition,Nothing}=nothing,
    DIR::AbstractString=@__DIR__,
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

function SolverSet(
    SET::AbstractSetup,
    PS::AbstractPhysicalSpace,
    GAS::AbstractProperty,
    IB::Union{AbstractCondition,Nothing}=nothing,
    DIR::AbstractString=@__DIR__,
)
    return SolverSet{typeof(SET),typeof(PS),Nothing,typeof(GAS),typeof(IB),typeof(DIR)}(
        SET,
        PS,
        PS,
        nothing,
        nothing,
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
    set = set_setup(; config...)

    # physical space
    ps = set_geometry(; config...)

    # velocity space
    vs = set_velocity(; config...)

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
        cp(file, string(outputFolder, "config.txt"))
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
