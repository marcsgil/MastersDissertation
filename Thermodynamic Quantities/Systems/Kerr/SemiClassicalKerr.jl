module SemiClassicalKerr

export sc_kerr


include("../../Utils.jl")
include("../../ThermodynamicSystems.jl")

using .Utils
using FastGaussQuadrature:gaussradau

#Auxiliary Functions

ϕ(u,τ,χ)= τ*(1+2*χ*u)

σ(u,τ,χ,ϕ) = u*sinh(ϕ) - τ*χ*u^2

sqrt_term(u,τ,χ,ϕ) = √abs(1+2*τ*χ*u*tanh(0.5*ϕ))

#Partition Function

integrand_Z(u,τ,χ,ϕ) = exp(-σ(u,τ,χ,ϕ) + 0.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) ) * sqrt_term(u,τ,χ,ϕ)

integrand_Z(u,τ,χ)= integrand_Z(u,τ,χ,ϕ(u,τ,χ))

rescaled_integrand_Z(u,(τ,χ,cutoff)) = 0.5*cutoff*integrand_Z(0.5*cutoff*(u+1),τ,χ)

function Z( τ, χ, cutoff_pars, nws)
    evaluate_integral( rescaled_integrand_Z, (τ,χ,find_cutoff(u->integrand_Z(u,τ,χ),cutoff_pars)), nws )
end

#Energy
function integrand_U(u,τ,χ,ϕ) 
    (           u*exp(-σ(u,τ,χ,ϕ) + 1.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^3 
         +  χ*u^2*exp(-σ(u,τ,χ,ϕ) + 2.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^5 ) * sqrt_term(u,τ,χ,ϕ)
 end
 
integrand_U(u,τ,χ) = integrand_U(u,τ,χ,ϕ(u,τ,χ))
 
rescaled_integrand_U(u,(τ,χ,cutoff)) = 0.5*cutoff*integrand_U(0.5*cutoff*(u+1),τ,χ)

function U(τ, χ, Z, cutoff_pars, nws)
    evaluate_integral( rescaled_integrand_U, (τ,χ,find_cutoff(u->integrand_U(u,τ,χ),cutoff_pars)), nws )/Z - 0.25χ
end

function U(τ, χ, cutoff_pars, nws)
    U(τ, χ, Z(τ, χ, cutoff_pars, nws), cutoff_pars, nws)
end

#Square of the energy
function integrand_U2(u,τ,χ,ϕ) 
    (             -2.5*χ*u*exp(-σ(u,τ,χ,ϕ) + 1.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^3 
         + (1-3.5*χ^2)*u^2*exp(-σ(u,τ,χ,ϕ) + 2.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^5 
                 + 2*χ*u^3*exp(-σ(u,τ,χ,ϕ) + 3.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^7 
                 + χ^2*u^4*exp(-σ(u,τ,χ,ϕ) + 4.5*ϕ) * ( 0.5*( 1+exp(-ϕ) ) )^9 ) * sqrt_term(u,τ,χ,ϕ)
end
 
integrand_U2(u,τ,χ) = integrand_U2(u,τ,χ,ϕ(u,τ,χ))
 
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


sc_kerr = ThermodynamicalSystem(( τ,χ )->  U( τ, χ, cutoff_pars, nws),
                                ( τ,χ )->  C( τ, χ, cutoff_pars, nws) )

end

using .SemiClassicalKerr
using BenchmarkTools
sc_kerr.U(.01,.1)
@benchmark sc_kerr.U(.01,.1)

using ThreadedIterables
θs = LinRange(.1,10,100)
@benchmark tmap(θ-> sc_kerr.U(θ,.1), θs)