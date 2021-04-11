using PyCall, Conda

# try importing meshio
# catch the installation
# 1) Julia built-in miniconda
# 2) global pip installer

cmd = `pip3 install meshio --user`
run(cmd)

#Conda.add_channel("conda-forge")
#Conda.add("meshio")

meshio = pyimport("meshio")
