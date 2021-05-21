using OrdinaryDiffEq, Plots
using KitBase

begin
    case = "homogeneous"
    maxTime = 3
    tlen = 16
    u0 = -5
    u1 = 5
    nu = 80
    nug = 0
    v0 = -5
    v1 = 5
    nv = 28
    nvg = 0
    w0 = -5
    w1 = 5
    nw = 28
    nwg = 0
    vMeshType = "rectangle"
    nm = 5
    knudsen = 1
    inK = 0
    alpha = 1.0
    omega = 0.5
    nh = 8
end