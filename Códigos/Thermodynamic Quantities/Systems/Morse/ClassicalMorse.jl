module ClassicalMorse

include("../../ThermodynamicSystems.jl")
include("../../Utils.jl")

using SpecialFunctions: dawson

export cl_morse

Z(τ,χ) = dawson(√(τ/4χ))/√τ

U,C = thermo_quantities_from_partition(Z)

cl_morse = ThermodynamicalSystem(U,C)

end 