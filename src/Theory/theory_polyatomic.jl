# ------------------------------------------------------------
# Rykov model
# ------------------------------------------------------------

"""
$(SIGNATURES)

Calculate dimensionless rotation number in Rykov model
"""
rykov_zr(Tₜ, T₀, Z₀) = Z₀ / (1.0 + π^1.5 / 2.0 * sqrt(T₀ / Tₜ) + (π + 0.25 * π^2) * T₀ / Tₜ)
