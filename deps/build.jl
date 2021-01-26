using PyCall

# try importing meshio
# if not existed, install it
try
    using Conda
    Conda.add_channel("conda-forge")
    Conda.add("meshio")
    meshio = pyimport("meshio")
catch
    cmd = `pip3 install meshio --user`
    run(cmd)
    meshio = pyimport("meshio")
end
