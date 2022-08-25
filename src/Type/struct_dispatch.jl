# ============================================================
# Structs for Multiple Dispatch
# ============================================================

"""
$(TYPEDEF)

Generic dispatcher

# Parameters
* `N`: type index
"""
struct Class{N} end


"""
$(TYPEDEF)

Velocity distribution marker for multiple dispatch

# Parameters
* `NF`: number of distribution functions
* `NV`: number of velocities
"""
struct VDF{NF,NV} end
