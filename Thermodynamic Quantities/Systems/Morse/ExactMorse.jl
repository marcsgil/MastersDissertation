module ExactMorse

export ex_morse

include("../../ThermodynamicSystems.jl")

Z(τ,χ) = partition_from_spectrum( (n,χ)-> (n+0.5)-χ*(n+0.5)^2, τ, χ, floor(Int64, 0.5*(1/χ-1)), 0.0 )

Z(1.0,1.0)

U,C = thermo_quantities_from_partition(Z)

ex_morse = ThermodynamicalSystem(U,C)

end