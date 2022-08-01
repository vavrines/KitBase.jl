# ============================================================
# Structs for Multiple Dispatch
# ============================================================

"""
$(TYPEDEF)

Velocity distribution marker for multiple dispatch

# Parameters
* `NF`: number of distribution functions
* `NV`: number of velocities
"""
struct VDF{NF,NV} end
