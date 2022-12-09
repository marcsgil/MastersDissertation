module EXKerr

include("../../build_thermo_systems.jl")

Z(θ,χ) = partition_from_spectrum( (n,χ)-> (n+0.5)+χ*(n+0.5)^2, θ, χ,  10^8, 10^-8 )

U = energy_from_partition(Z)
C = heat_from_partition(Z)

end