macro nametuple(x...)
    ex = [:($(esc(z)) = $(esc(z))) for z in x]
    return :(return ($(ex...),))
end
