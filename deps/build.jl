using PyCall, Conda

# try importing
# catch the installation
# 1) Julia built-in miniconda
# 2) global pip installer

#try
#    using Conda
#    Conda.add_channel("conda-forge")
#    Conda.add("meshio")
#    Conda.add("scipy")
#    pyimport("meshio")
#    pyimport("scipy")
#catch
    cmd = `pip3 install meshio scipy --user`
    run(cmd)
#    pyimport("meshio")
#    pyimport("scipy")
#end
