const itp = PyNULL()

function __init__()
    #=np = nworkers()
    nt = Threads.nthreads()

    show_worker(np) = begin
        if np == 1
            "$np worker"
        else
            "$np workers"
        end
    end
    show_thread(nt) = begin
        if nt == 1
            "$nt thread"
        else
            "$nt threads"
        end
    end

    if has_cuda()
        @info "Kinetic will run with $(show_worker(np)), $(show_thread(nt)) and CUDA"
        for (i, dev) in enumerate(CUDA.devices())
            @info "$i: $(CUDA.name(dev))"
        end
        #@info "Scalar operation is disabled in CUDA"
        CUDA.allowscalar(false)
    else
        @info "Kinetic will run with $(show_worker(np)) and $(show_thread(nt))"
    end=#

    copy!(itp, pyimport("scipy.interpolate"))
end

"""
$(SIGNATURES)

RHS-ODE of Boltzmann equation with non-uniform velocity
"""
function boltzmann_nuode!(df, f::AA{T,3}, p, t) where {T}
    Kn, M, phi, psi, phipsi, u, v, w, vnu, u1, v1, w1, vuni = p

    nu = length(u)
    nv = length(v)
    nw = length(w)
    nu1 = length(u1)
    nv1 = length(v1)
    nw1 = length(w1)

    curve = itp.RegularGridInterpolator((u, v, w), f)
    _f = reshape(curve(vuni), nu1, nv1, nw1)
    _df = boltzmann_fft(_f, Kn, M, phi, psi, phipsi)

    curve1 = itp.RegularGridInterpolator((u1, v1, w1), _df)
    df .= reshape(curve1(vnu), nu, nv, nw)
end

boltzmann_nuode!(
    zeros(16, 16, 16),
    rand(16, 16, 16),
    (5.0, 5, phi, psi, phipsi, u, v, w, vnu, uuni1d, vuni1d, wuni1d, vuni),
    0.0,
)

"""
$(SIGNATURES)

Maxwell quadrature

## Arguments
* `N`: quadrature order (MUST less than 33)
"""
function maxwell_quadrature(N::Integer, C = 1)
    @assert N <= 33

    py"""
    import numpy as np
    from numpy import linalg as LA

    def dvGH(N2,C):
        N = N2//2

        a = np.zeros(N)
        b = np.zeros(N)
        a[0] = 1.0 / np.sqrt(np.pi)
        a[1] = 2.0 / np.sqrt(np.pi) / (np.pi - 2.0)
        b[1] = a[0] / (a[0] + a[1]) / 2.0

        for i in range(2,N):
            b[i] = (i - 1) + 0.5 - b[i-1] - a[i-1]**2
            a[i] = (i**2 / 4.0 / b[i] - b[i-1] - 0.5) / a[i-1] - a[i-1]

        J = np.diag(a) + np.diag(np.sqrt(b[1:N]), 1) \
        + np.diag(np.sqrt(b[1:N]), -1)

        v, V = LA.eig(J)

        w = V[0, :] * V[0, :] * np.sqrt(np.pi) / 2.0

        vw = np.transpose(np.vstack((v, w)))
        vw = vw[vw[:, 0].argsort()]
        v = vw[:, 0]
        w = vw[:, 1]

        Xis = np.hstack((-np.flipud(v), v))
        weights = np.hstack((np.flipud(w), w))
        weights = weights * np.exp(Xis**2) * C
        Xis = Xis*C
        return (Xis, weights)
    """

    p, w = py"dvGH"(N, C)
end
