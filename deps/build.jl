using PyCall, Conda
import Pkg

# try importing meshio
# catch the installation
# 1) Julia built-in miniconda
# 2) global pip installer

JULIA_KIT_DIR = homedir()

#ENV["PYTHON"] = "/usr/bin/python"
ENV["PYTHON"] = "$(JULIA_KIT_DIR)/.julia/conda/3/bin/python"
Pkg.build("PyCall")

PIP = "$(JULIA_KIT_DIR)/.julia/conda/3/bin/pip"
run(`$PIP install meshio --user`)

#Conda.add_channel("conda-forge")
#Conda.add("meshio")
#Conda.pip_interop(true)
#Conda.pip("install", "meshio")

#@info "installing meshio"
#cmd = `pip install meshio --user`
#cmd = `./usr/bin/pip install meshio --user`
#run(cmd)

meshio = pyimport("meshio")
