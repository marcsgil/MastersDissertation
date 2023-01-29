module ExactPöschlTeller

export ex_pt

include("../../ThermodynamicSystems.jl")

Z(τ,χ) = partition_from_spectrum( (n,χ)-> (n+0.5)-χ*(n+0.5)^2-0.25*χ, τ, χ, floor(Int64, 0.5*(1/χ-1)), 0.0 )

U,C = thermo_quantities_from_partition(Z)

ex_pt = ThermodynamicalSystem(U,C)

end