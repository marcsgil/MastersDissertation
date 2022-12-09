module DampenedMetaplecticKerr

export dm_kerr

include("../../Utils.jl")
include("../../ThermodynamicSystems.jl")

using .Utils
using FastGaussQuadrature:gaussradau

#Auxiliary Stuff
ϕ(u,τ,χ) = 0.5*τ*√( (1+2*χ*u)*(1+6*χ*u) )

f(ϕ) = 1-tanh(ϕ)/ϕ

σ(u,χ,ϕ) = u+χ*u^2-f(ϕ)*u*(1+2*χ*u)^2/(1+6*χ*u)

#Partition Function

integrand_Z(u,τ,χ,ϕ) = sech(ϕ)*exp( -τ*σ(u,χ,ϕ) )

integrand_Z(u,τ,χ) = integrand_Z(u,τ,χ,ϕ(u,τ,χ))

rescaled_integrand_Z(u,(τ,χ,cutoff)) = 0.5*cutoff*integrand_Z(0.5*cutoff*(u+1),τ,χ)

function Z( τ, χ, cutoff_pars, nws)
    evaluate_integral( rescaled_integrand_Z, (τ,χ,find_cutoff(u->integrand_Z(u,τ,χ),cutoff_pars)), nws )
end

#Energy

integrand_U(u,τ,χ) = (u+χ*u^2)*integrand_Z(u,τ,χ)

rescaled_integrand_U(u,(τ,χ,cutoff)) = 0.5*cutoff*integrand_U(0.5*cutoff*(u+1),τ,χ)

function U(τ, χ, Z, cutoff_pars, nws)
    evaluate_integral( rescaled_integrand_U, (τ,χ,find_cutoff(u->integrand_U(u,τ,χ),cutoff_pars)), nws )/Z
end

function U(τ, χ, cutoff_pars, nws)
    U(τ, χ, Z(τ, χ, cutoff_pars, nws), cutoff_pars, nws)
end

#Square of the Energy
integrand_U2(u,τ,χ) = ( -2.5*χ*u + (1-3.5*χ^2)*u^2 + 2*χ*u^3 + χ^2*u^4)*integrand_Z(u,τ,χ)

rescaled_integrand_U2(u,(τ,χ,cutoff)) = 0.5*cutoff*integrand_U2(0.5*cutoff*(u+1),τ,χ)

function U2(τ, χ, Z, cutoff_pars, nws)
    evaluate_integral( rescaled_integrand_U2, (τ,χ,find_cutoff(u->integrand_U2(u,τ,χ),cutoff_pars)), nws )/Z + 5*χ^2/16 - 0.25
end

#Heat Capacity
C(τ, χ, Z, cutoff_pars, nws) = τ^2*( U2(τ, χ, Z, cutoff_pars, nws) - U(τ, χ, Z, cutoff_pars, nws)^2 )
C(τ, χ, cutoff_pars, nws) = C(τ, χ, Z( τ, χ, cutoff_pars, nws), cutoff_pars, nws)

#Creating the system

const nws = gaussradau(20)

#Cutoff Stuff
#(tol,0,max_cutoff,ϵ)
const cutoff_pars = (10^-5,0.0,5000.0,10^-8)


dm_kerr = ThermodynamicalSystem(( τ,χ )->  U( τ, χ, cutoff_pars, nws),
                                ( τ,χ )->  C( τ, χ, cutoff_pars, nws) )

end