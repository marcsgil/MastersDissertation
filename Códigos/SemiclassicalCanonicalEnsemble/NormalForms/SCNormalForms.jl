using DifferentialEquations,StaticArrays,Polynomials

function Z_integrand(J,θ,F::Polynomial,ω::Polynomial,ω′::Polynomial)
    #Calculates the integrand of the partition function
    ϕ = θ*ω(J)
    half_ϕ = ϕ/2
    Sᴱ = (ϕ-sinh(ϕ))*J - θ*F(J)
    modified_cosh = ( 1+exp(-ϕ) )/2
    sqrt_term = √abs((1+J*ω′(J)*θ*tanh(half_ϕ)))
    sqrt_term*modified_cosh*exp(Sᴱ + half_ϕ)
end

function A_integrand(J,θ,F::Polynomial,ω::Polynomial,ω′::Polynomial,G::Polynomial)
    #Calculates the integrand of the expectation value of an operator A = G(J)
    ϕ = θ*ω(J)
    half_ϕ = ϕ/2
    Sᴱ = (ϕ-sinh(ϕ))*J - θ*F(J)
    modified_cosh = ( 1+exp(-ϕ) )/2
    sqrt_term = √abs((1+J*ω′(J)*θ*tanh(half_ϕ)))
    sqrt_term*sum( pair-> pair[2]*J^(pair[1])*modified_cosh^(2*pair[1]+1)*exp(Sᴱ + (2*pair[1]+1)*half_ϕ), pairs(G))
end

function ode_step(y,(θ,F,ω,ω′,G),J)
    dy1 = Z_integrand(J,θ,F,ω,ω′)
    dy2 = A_integrand(J,θ,F,ω,ω′,G)
    SA[dy1,dy2]
end

function expectation_value( θ,F::Polynomial,G::Polynomial,J_max=200.)
    ω = derivative(F)
    ω′ = derivative(ω)
    prob = ODEProblem(ode_step,SA[zero(J_max),zero(J_max)],(0,J_max),(θ,F,ω,ω′,G))
    sol = solve(prob,save_everystep=false,save_start=false)
    sol[2,end]/sol[1,end]
end

function expectation_value( θs::AbstractArray,F::Polynomial,G::Polynomial,J_max=200.)
    ω = derivative(F)
    ω′ = derivative(ω)
    prob = ODEProblem(ode_step,zeros(SVector{2,typeof(J_max)}),(0,J_max),(θs[1],F,ω,ω′,G))

    function prob_func(prob,i,repeat)
        remake(prob,p=(θs[i],F,ω,ω′,G))
    end

    output_func(sol,i) = (sol[2,end]/sol[1,end],false)

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)
    solve(ensemble_prob,trajectories=length(θs))
end