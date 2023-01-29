using DifferentialEquations,StaticArrays,Polynomials

fr(z) = ifelse(iszero(z),one(typeof(z)),tanh(z)/z)
fi(z) = ifelse(iszero(z),one(typeof(z)),tan(z)/z)

f(α) = real(ifelse(isreal(α),fr(α),fi(imag(α))))
g(α) = real((1-f(α))/α^2)


function Z_integrand(J,θ,F,ω::Polynomial,Ω²::Polynomial)
    #Calculates the integrand of the partition function
    α = (θ/2)*√Complex(Ω²(J))
    Real(sech(α))*exp(θ*( J*ω(J)^3*g(α)*θ^2/4- F(J) ))
end

function ode_step(y,(θ,F,ω,Ω²,G),J)
    Z_integrand(J,θ,F,ω,Ω²)*SA[1,G(J)]
end

function expectation_value( θs::AbstractArray,F::Polynomial,G::Polynomial,J_max=200.)
    ω = derivative(F)
    ω′ = derivative(ω)
    pol_J = Polynomial([0,1])
    Ω² = ω*(ω+2*pol_J*ω′)
    prob = ODEProblem(ode_step,zeros(SVector{2,typeof(J_max)}),(0,J_max),(θs[1],F,ω,Ω²,G))

    function prob_func(prob,i,repeat)
        remake(prob,p=(θs[i],F,ω,Ω²,G))
    end

    output_func(sol,i) = (sol[2,end]/sol[1,end],false)

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)
    solve(ensemble_prob,trajectories=length(θs))
end