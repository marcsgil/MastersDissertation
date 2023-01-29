module ClassicalKerr

include("../../ThermodynamicSystems.jl")

using SpecialFunctions

export cl_kerr

U(τ,χ) =  0.25*( 2/τ -1/χ + 2*exp(-τ/(4*χ))/( √(π*χ*τ)*erfc(0.5*√(τ/χ)) ) )

cl_kerr = ThermodynamicalSystem(U)

end 