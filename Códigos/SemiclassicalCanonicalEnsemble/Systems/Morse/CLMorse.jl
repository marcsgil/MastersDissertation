module CLMorse

include("../../build_thermo_systems.jl")

using SpecialFunctions: dawson

Z(θ,χ) = dawson(√(θ/4χ))/√θ

U = energy_from_partition(Z)
C = heat_from_partition(Z)
end 