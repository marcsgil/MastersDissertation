module CLKerr

include("../../build_thermo_systems.jl")
using SpecialFunctions

function Z(θ,χ)
    par = θ/(4χ)
    √(π/(θ*χ))*exp(par)*erfc(√par)/2
end

U = energy_from_partition(Z)
C = heat_from_partition(Z)

end 