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

Dimension dispatcher

# Parameters
* `N`: dimension index
"""
struct Dimension{N} end


"""
$(TYPEDEF)

Velocity distribution function dispatcher

# Parameters
* `NF`: number of distribution functions
* `NV`: number of velocities
"""
struct VDF{NF,NV} end
