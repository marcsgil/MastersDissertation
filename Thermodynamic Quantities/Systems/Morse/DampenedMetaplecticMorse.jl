module DampenedMetaplecticMorse

export dm_morse

__revise_mode__ = :eval

include("../../Utils.jl")
include("../../ThermodynamicSystems.jl")

using .Utils
using FastGaussQuadrature, IfElse
using SpecialFunctions:erf

#Auxiliar Functions
ϕ(y,τ) = 0.5*τ*√complex((1-2y)*(1-y))

f(ϕ) = IfElse.ifelse(iszero(ϕ),1.0,tanh(ϕ)/ϕ)

g(α) = IfElse.ifelse(iszero(α),2/√π,erf(α)/α)

h(ϕ) = IfElse.ifelse(iszero(ϕ),1/3,(1-f(ϕ))/ϕ^2)

σ(y,τ,χ,ϕ) = τ*y^2*( 1-0.25*τ^2*(1-y)^2*h(ϕ) )/(4χ)

α(y,τ,χ,ϕ) = complex(τ*f(ϕ)*(1-y^2)/(4χ))

#Partition Function

integrand_Z(y,τ,χ,ϕ) = sech(ϕ)*exp(-σ(y,τ,χ,ϕ))*g(√α(y,τ,χ,ϕ))

integrand_Z(y,(τ,χ)) = real(integrand_Z(y,τ,χ,ϕ(y,τ)))

nws = gausschebyshev(800,3)

Z(τ,χ) = evaluate_complex_integral(integrand_Z,(τ,χ),nws)

#Energy

function integrand_U(y,τ,χ,ϕ,f,α)
    sech(ϕ)*exp(-σ(y,τ,χ,ϕ))*( g(√α)*(1/(2*τ*f) +y^2/(4χ)) -exp(-α)/(τ*√π*f) )
end

integrand_U(y,τ,χ,ϕ) = integrand_U(y,τ,χ,ϕ,f(ϕ),α(y,τ,χ,ϕ))

integrand_U(y,(τ,χ)) = real(integrand_U(y,τ,χ,ϕ(y,τ)))

U(τ,χ) = evaluate_complex_integral(integrand_U,(τ,χ),nws)/Z(τ,χ)

#Heat Capacity
C(τ,χ) = 0

dm_morse = ThermodynamicalSystem(U,C)

#χs for some diatomic molecules
χH2 = 121.34/4401.21
χN2 = 14.32/2358.56
χO2 = 11.98/1580.19
χCs2=0.08/42.02

end