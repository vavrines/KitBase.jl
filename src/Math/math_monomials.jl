degree_size(currDegree, spatialDim) =
    factorial(currDegree + spatialDim - 1) ÷
    (factorial(currDegree) * factorial(spatialDim - 1))


function basis_size(LMaxDegree, spatialDim)
    basisLen = 0
    for idx_degree = 0:LMaxDegree
        basisLen += degree_size(idx_degree, spatialDim)
    end
    return basisLen
end


function monomial_basis(pointX, pointY, pointZ, polyDegree::T) where {T<:Integer}
    idx_vector = 1
    spatialDim = 3
    basisLen = basis_size(polyDegree, spatialDim)
    basisAtPt = ones(1, basisLen)
    for idx_degree = 0:polyDegree
        for a = 0:idx_degree
            for b = 0:(idx_degree-a)
                c = idx_degree - a - b
                basisAtPt[idx_vector] = pointX^a * pointY^b * pointZ^c
                idx_vector = idx_vector + 1
            end
        end
    end

    return basisAtPt
end


function monomial_basis(pointX, polyDegree::T) where {T<:Integer}
    idx_vector = 1
    spatialDim = 1
    basisLen = basis_size(polyDegree, spatialDim)
    basisAtPt = ones(1, basisLen)
    for a = 0:polyDegree
        basisAtPt[idx_vector] = pointX^a
        idx_vector = idx_vector + 1
    end

    return basisAtPt
end


function eval_sphermonomial(quadpts::AV, polyDegree::Integer)
    monomialBasis = zeros(basis_size(polyDegree, 1), size(quadpts, 1))

    for idx_quad = 1:length(quadpts)
        monomialBasis[:, idx_quad] = monomial_basis(quadpts[idx_quad, 1], polyDegree)
    end

    return monomialBasis
end


function eval_sphermonomial(quadpts::Matrix, polyDegree::Integer)
    monomialBasis = zeros(basis_size(polyDegree, 3), size(quadpts, 1))

    for idx_quad = 1:(size(quadpts)[1])
        monomialBasis[:, idx_quad] = monomial_basis(
            quadpts[idx_quad, 1],
            quadpts[idx_quad, 2],
            quadpts[idx_quad, 3],
            polyDegree,
        )
    end

    return monomialBasis
end
