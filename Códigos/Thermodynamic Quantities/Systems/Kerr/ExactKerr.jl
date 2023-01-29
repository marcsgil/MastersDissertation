module ExactKerr

include("../../ThermodynamicSystems.jl")

export  ex_kerr

U,C = thermo_quantities_from_spectrum( (n, χ) -> (n+0.5)*( 1+χ*(n+0.5) ), 10^8,10^-8 )

ex_kerr = ThermodynamicalSystem(U,C)

end