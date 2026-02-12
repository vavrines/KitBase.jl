"""
$(SIGNATURES)

Workspace for 0F 1F 2F step!
"""
struct StepWorkspace{T1<:AV, T2, T3} <: AbstractWorkspace
    w_old::T1
    
    MH::T2
    SH::T2
    MB::T3
    SB::T3
end

function StepWorkspace(ctr::ControlVolume)
    return StepWorkspace(
        similar(ctr.w),
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function StepWorkspace(ctr::ControlVolume1F)
    return StepWorkspace(
        similar(ctr.w),
        similar(ctr.f),
        similar(ctr.f),
        nothing,
        nothing,
    )
end

function StepWorkspace(ctr::ControlVolume2F)
    return StepWorkspace(
        similar(ctr.w),
        similar(ctr.h),
        similar(ctr.h),
        similar(ctr.b),
        similar(ctr.b),
    )
end

"""
$(SIGNATURES)

Workspace for 3F Rykov step!
"""
struct StepWorkspace3F{T1<:AV, T2, T3, T4} <: AbstractWorkspace
    w_old::T1
    
    MHT::T2
    MHR::T2
    MH::T2
    SHT::T2
    SHR::T2
    MBT::T3
    MBR::T3
    MB::T3
    SBT::T3
    SBR::T3
    MRT::T4
    MRR::T4
    MR::T4
    SRT::T4
    SRR::T4
end

function StepWorkspace3F(ctr::ControlVolume3F)
    return StepWorkspace3F(
        similar(ctr.w),
        similar(ctr.h0),
        similar(ctr.h0),
        similar(ctr.h0),
        similar(ctr.h0),
        similar(ctr.h0),
        similar(ctr.h1),
        similar(ctr.h1),
        similar(ctr.h1),
        similar(ctr.h1),
        similar(ctr.h1),
        similar(ctr.h2),
        similar(ctr.h2),
        similar(ctr.h2),
        similar(ctr.h2),
        similar(ctr.h2),
    )
end


"""
$(SIGNATURES)

Workspace for 0F 1F 2F 3F mixture step!
"""
struct StepWorkspaceM{T1<:AM, T2<:AV, T3, T4, T5, T6} <: AbstractWorkspace
    w_old::T1
    prim_old::T1

    mw::T1
    mprim::T1

    tau::T2
    
    H0::T3
    H1::T4
    H2::T5
    H3::T6
end

function StepWorkspaceM(ctr::ControlVolume)
    return StepWorkspaceM(
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w[1, :]),
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function StepWorkspaceM(ctr::ControlVolume1F)
    return StepWorkspaceM(
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w[1, :]),
        similar(ctr.f),
        nothing,
        nothing,
        nothing,
    )
end

function StepWorkspaceM(ctr::ControlVolume2F)
    return StepWorkspaceM(
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w[1, :]),
        similar(ctr.h),
        similar(ctr.b),
        nothing,
        nothing,
    )
end

function StepWorkspaceM(ctr::ControlVolume3F)
    return StepWorkspaceM(
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w[1, :]),
        similar(ctr.h0),
        similar(ctr.h1),
        similar(ctr.h2),
        nothing,
    )
end

function StepWorkspaceM(ctr::ControlVolume4F)
    return StepWorkspaceM(
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w),
        similar(ctr.w[1, :]),
        similar(ctr.h0),
        similar(ctr.h1),
        similar(ctr.h2),
        similar(ctr.h3),
    )
end