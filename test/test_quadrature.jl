# spherical quadrature
KitBase.legendre_quadrature(6)
KitBase.octa_quadrature(8)

KitBase.VSpace1D() |> show
KitBase.VSpace2D() |> show
KitBase.VSpace3D() |> show
KitBase.MVSpace1D() |> show
KitBase.MVSpace2D() |> show
KitBase.MVSpace3D() |> show

KitBase.VSpace1D(-5, 5, 16, "newton")
KitBase.VSpace2D(-5, 5, 16, -5, 5, 16, "newton")
KitBase.VSpace2D(-5, 5, 16, -5, 5, 16, "algebra")
KitBase.VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16, "newton")
KitBase.MVSpace1D(-5, 5, -5, 5, 16, "newton")
KitBase.MVSpace2D(-5, 5, -5, 5, 16, -5, 5, -5, 5, 16, "newton")
KitBase.MVSpace3D(-5, 5, -5, 5, 8, -5, 5, -5, 5, 8, -5, 5, -5, 5, 8, "newton")
KitBase.UnstructVSpace(-1, 1, 16, rand(16), ones(16))
