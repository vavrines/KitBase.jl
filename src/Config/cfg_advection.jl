function ib_advection(args...)
    fw = (x, p...) -> sin(2π * x)
    return fw, nothing, NamedTuple()
end
