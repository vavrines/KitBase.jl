"""
$(SIGNATURES)

Calculate electromagnetic coeffcients in hyperbolic Maxwell's equations
"""
function em_coefficients(prim::AM, E::AV, B::AV, mr, lD, rL, dt)

    if eltype(prim) <: Int
        A = zeros(9, 9)
        b = zeros(9)
    else
        A = similar(prim, 9, 9)
        A .= 0.0
        b = similar(prim, 9)
        b .= 0.0
    end

    A[1, 1] = -1.0 / (2.0 * rL)
    A[2, 2] = -1.0 / (2.0 * rL)
    A[3, 3] = -1.0 / (2.0 * rL)
    A[4, 1] = mr / (2.0 * rL)
    A[5, 2] = mr / (2.0 * rL)
    A[6, 3] = mr / (2.0 * rL)
    A[7, 1] = 1.0 / (dt)
    A[8, 2] = 1.0 / (dt)
    A[9, 3] = 1.0 / (dt)

    A[1, 4] = 1.0 / (dt)
    A[1, 5] = -B[3] / (2.0 * rL)
    A[1, 6] = B[2] / (2.0 * rL)
    A[2, 4] = B[3] / (2.0 * rL)
    A[2, 5] = 1.0 / (dt)
    A[2, 6] = -B[1] / (2.0 * rL)
    A[3, 4] = -B[2] / (2.0 * rL)
    A[3, 5] = B[1] / (2.0 * rL)
    A[3, 6] = 1.0 / (dt)

    A[4, 7] = 1.0 / (dt)
    A[4, 8] = mr * B[3] / (2.0 * rL)
    A[4, 9] = -mr * B[2] / (2.0 * rL)
    A[5, 7] = -mr * B[3] / (2.0 * rL)
    A[5, 8] = 1.0 / (dt)
    A[5, 9] = mr * B[1] / (2.0 * rL)
    A[6, 7] = mr * B[2] / (2.0 * rL)
    A[6, 8] = -mr * B[1] / (2.0 * rL)
    A[6, 9] = 1.0 / (dt)

    A[7, 4] = prim[1, 1] / (2.0 * rL * lD^2)
    A[8, 5] = prim[1, 1] / (2.0 * rL * lD^2)
    A[9, 6] = prim[1, 1] / (2.0 * rL * lD^2)
    A[7, 7] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[8, 8] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[9, 9] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)

    b[1] =
        prim[2, 1] / (dt) + E[1] / (2.0 * rL) - B[2] * prim[4, 1] / (2.0 * rL) +
        B[3] * prim[3, 1] / (2.0 * rL)
    b[2] =
        prim[3, 1] / (dt) + E[2] / (2.0 * rL) - B[3] * prim[2, 1] / (2.0 * rL) +
        B[1] * prim[4, 1] / (2.0 * rL)
    b[3] =
        prim[4, 1] / (dt) + E[3] / (2.0 * rL) - B[1] * prim[3, 1] / (2.0 * rL) +
        B[2] * prim[2, 1] / (2.0 * rL)
    b[4] =
        prim[2, 2] / (dt) - mr * E[1] / (2.0 * rL) + mr * B[2] * prim[4, 2] / (2.0 * rL) -
        mr * B[3] * prim[3, 2] / (2.0 * rL)
    b[5] =
        prim[3, 2] / (dt) - mr * E[2] / (2.0 * rL) + mr * B[3] * prim[2, 2] / (2.0 * rL) -
        mr * B[1] * prim[4, 2] / (2.0 * rL)
    b[6] =
        prim[4, 2] / (dt) - mr * E[3] / (2.0 * rL) + mr * B[1] * prim[3, 2] / (2.0 * rL) -
        mr * B[2] * prim[2, 2] / (2.0 * rL)
    b[7] =
        E[1] / (dt) - prim[1, 1] * prim[2, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[2, 2] * mr / (2.0 * rL * lD^2)
    b[8] =
        E[2] / (dt) - prim[1, 1] * prim[3, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[3, 2] * mr / (2.0 * rL * lD^2)
    b[9] =
        E[3] / (dt) - prim[1, 1] * prim[4, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[4, 2] * mr / (2.0 * rL * lD^2)

    return A, b

end
