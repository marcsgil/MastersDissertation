module ClassicalPöschlTeller

include("../../ThermodynamicSystems.jl")
include("../../Utils.jl")

using SpecialFunctions, FastGaussQuadrature, IfElse, .Utils

export cl_pt

g(α) = IfElse.ifelse(iszero(α),2/√π,erf(α)/α)

integrand(Q,(τ,χ)) = exp(-0.25*τ*(1/χ-χ)Q^2)*g(√( 0.25*τ*(1/χ-χ)*(1-Q^2) ))

const nws = gausschebyshev(100)

Z(τ,χ) = evaluate_integral(integrand, (τ,χ), nws)

U,C = thermo_quantities_from_partition(Z)

cl_pt = ThermodynamicalSystem(U,C)

end 