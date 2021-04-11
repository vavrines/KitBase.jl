using PyCall, Conda
import Pkg

# try importing meshio
# catch the installation
# 1) Julia built-in miniconda
# 2) global pip installer

ENV["PYTHON"] = "/usr/bin/python3"
Pkg.build("PyCall")

#Conda.add_channel("conda-forge")
#Conda.add("meshio")
#Conda.pip_interop(true)
#Conda.pip("install", "meshio")

@info "installing meshio"
cmd = `pip3 install meshio --user`
run(cmd)

meshio = pyimport("meshio")
