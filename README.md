# KitBase.jl

[![version](https://juliahub.com/docs/KitBase/version.svg)](https://juliahub.com/ui/Packages/KitBase/YOFTS)
![CI](https://github.com/vavrines/KitBase.jl/workflows/CI/badge.svg)
[![](https://img.shields.io/badge/docs-stable-green.svg)](https://xiaotianbai.com/Kinetic.jl/stable/)
[![codecov](https://codecov.io/gh/vavrines/KitBase.jl/branch/main/graph/badge.svg?token=vGgQhyGJ6L)](https://codecov.io/gh/vavrines/KitBase.jl)
[![deps](https://juliahub.com/docs/KitBase/deps.svg)](https://juliahub.com/ui/Packages/KitBase/YOFTS?t=2)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/vavrines/KitBase.jl/v0.8.0.svg?style=social&logo=github)](https://github.com/vavrines/KitBase.jl)

This lightweight module provides basic physical formulations and numerical methods for [Kinetic.jl](https://github.com/vavrines/Kinetic.jl) ecosystem.
The finite volume method (FVM) is employed to perform 1-3 dimensional numerical simulations on CPUs and GPUs.
A rich set of numerical fluxes and source terms are implemented for advection-diffusion-type governing equations.
A partial list of current supported models and equations include:
- Boltzmann equation
- radiative transfer equation
- Fokker-Planck-Landau equation
- direct simulation Monte Carlo
- advection-diffusion equation
- Burgers equation
- Euler equations
- Navier-Stokes equations
- Magnetohydrodynamical equations
- Maxwell's equations

For the detailed information on the implementation and usage of the package, you may
[check the documentation](https://xiaotianbai.com/Kinetic.jl/dev/).
