module EXMorse

include("../../build_thermo_systems.jl")

Z(τ,χ) = partition_from_spectrum( (n,χ)-> (n+0.5)-χ*(n+0.5)^2, τ, χ, floor(Int64, 0.5*(1/χ-1)), 0.0 )

U = energy_from_partition(Z)
C = heat_from_partition(Z)
end