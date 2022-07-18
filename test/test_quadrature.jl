# spherical quadrature
KitBase.legendre_quadrature(6)
KitBase.octa_quadrature(8)

KitBase.vs1D() |> show
KitBase.vs2D() |> show
KitBase.vs3D() |> show
KitBase.MVSpace1D() |> show
KitBase.MVSpace2D() |> show
KitBase.MVSpace3D() |> show

KitBase.vs1D(-5, 5, 16, "newton")
KitBase.vs1D(-5, 5, 16, "algebra")
KitBase.vs2D(-5, 5, 16, -5, 5, 16, "newton")
KitBase.vs2D(-5, 5, 16, -5, 5, 16, "algebra")
KitBase.vs2D(-5, 5, 16, -5, 5, 16, "maxwell")
KitBase.vs3D(-5, 5, 16, -5, 5, 16, -5, 5, 16, "newton")
KitBase.vs3D(-5, 5, 16, -5, 5, 16, -5, 5, 16, "algebra")
KitBase.MVSpace1D(-5, 5, -5, 5, 16, "newton")
KitBase.MVSpace2D(-5, 5, -5, 5, 16, -5, 5, -5, 5, 16, "newton")
KitBase.MVSpace3D(-5, 5, -5, 5, 8, -5, 5, -5, 5, 8, -5, 5, -5, 5, 8, "newton")
KitBase.UnstructVSpace(-1, 1, 16, rand(16), ones(16))
