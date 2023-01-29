module DMPöschlTeller

export dm_pt

__revise_mode__ = :eval

include("../../Utils.jl")
include("../../ThermodynamicSystems.jl")

using .Utils
using FastGaussQuadrature, IfElse
using SpecialFunctions:erf

#Auxiliar Functions
ϕ(Q,τ,χ) = 0.5*τ*(1-χ^2)*√complex((1-Q^2)*(1-3Q^2))

f(ϕ) = IfElse.ifelse(iszero(ϕ),zero(ϕ),tanh(ϕ)/ϕ)

g(α) = IfElse.ifelse(iszero(α),2/√π,erf(α)/α)

h(ϕ) = IfElse.ifelse(iszero(ϕ),1/3,(1-f(ϕ))/ϕ^2)

σ(Q,τ,χ,ϕ) = 0.25*τ*(1/χ-χ)*Q^2*( 1-0.25*τ^2*(1-χ^2)*(1-Q^2)^2*h(ϕ) )

α(Q,τ,χ,ϕ) = complex(0.25*τ*f(ϕ)*(1/χ-χ)*(1-Q^2))

#Partition Function

integrand_Z(Q,τ,χ,ϕ) = sech(ϕ)*exp(-σ(Q,τ,χ,ϕ))*g(√α(Q,τ,χ,ϕ))

integrand_Z(Q,(τ,χ)) = real(integrand_Z(Q,τ,χ,ϕ(Q,τ,χ)))

nws = gausschebyshev(100)

Z(τ,χ) = evaluate_complex_integral(integrand_Z,(τ,χ),nws)

#Energy

function integrand_U(Q,τ,χ,ϕ,f,α)
    sech(ϕ)*exp(-σ(Q,τ,χ,ϕ))*( g(√α)*(1/(2*τ*f) + 0.25*Q^2*(1/χ-χ)) -exp(-α)/(τ*√π*f) )
end

integrand_U(Q,τ,χ,ϕ) = integrand_U(Q,τ,χ,ϕ,f(ϕ),α(Q,τ,χ,ϕ))

integrand_U(Q,(τ,χ)) = real(integrand_U(Q,τ,χ,ϕ(Q,τ,χ)))

U(τ,χ) = evaluate_complex_integral(integrand_U,(τ,χ),nws)/Z(τ,χ)

#Heat Capacity
C(τ,χ) = 0

dm_pt = ThermodynamicalSystem(U,C)

#χs for some diatomic molecules
χH2 = 121.34/4401.21
χN2 = 14.32/2358.56
χO2 = 11.98/1580.19
χCs2=0.08/42.02

end