# spherical quadrature
KB.legendre_quadrature(6)
KB.octa_quadrature(8)

KB.VSpace1D() |> show
KB.VSpace2D() |> show
KB.VSpace3D() |> show
KB.MVSpace1D() |> show
KB.MVSpace2D() |> show
KB.MVSpace3D() |> show

KB.VSpace1D(-5, 5, 16; type="newton")
KB.VSpace1D(-5, 5, 16; type="algebra")
KB.VSpace2D(-5, 5, 16, -5, 5, 16; type="newton")
KB.VSpace2D(-5, 5, 16, -5, 5, 16; type="algebra")
KB.VSpace2D(-5, 5, 16, -5, 5, 16; type="maxwell")
KB.VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16; type="newton")
KB.VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16; type="algebra")
KB.MVSpace1D(-5, 5, -5, 5, 16; type="newton")
KB.MVSpace2D(-5, 5, -5, 5, 16, -5, 5, -5, 5, 16; type="newton")
KB.MVSpace3D(-5, 5, -5, 5, 8, -5, 5, -5, 5, 8, -5, 5, -5, 5, 8; type="newton")
KB.UnstructVSpace(-1, 1, 16, rand(16), ones(16))

mesh_quadrature(-5, -5, 8; type="rectangle")
mesh_quadrature(-5, -5, 9; type="newton")
mesh_quadrature(-5, -5, 8; type="algebra")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type="rectangle")
mesh_quadrature(-5, -5, 9, -5, -5, 9; type="newton")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type="algebra")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type="maxwell")
mesh_quadrature(-5, -5, 8, -5, -5, 8, -5, -5, 8; type="rectangle")
mesh_quadrature(-5, -5, 9, -5, -5, 9, -5, -5, 9; type="newton")
mesh_quadrature(-5, -5, 8, -5, -5, 8, -5, -5, 8; type="algebra")
