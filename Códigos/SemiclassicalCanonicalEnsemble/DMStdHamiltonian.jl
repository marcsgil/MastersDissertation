using DifferentialEquations,StaticArrays,SpecialFunctions
using ForwardDiff:derivative

fr(z) = ifelse(iszero(z),one(typeof(z)),tanh(z)/z)
fi(z) = ifelse(iszero(z),one(typeof(z)),tan(z)/z)

f(α) = real(ifelse(isreal(α),fr(α),fi(imag(α))))
g(α) = real((1-f(α))/α^2)

function I(θ,f,modified_V,n::Int)
    @assert n ≥ 0 "n must be greater than 0!"
    if n == 0
        erf(√(θ*f*modified_V))
    else
        ((n-1/2)*I(θ,f,modified_V,n-1) - modified_V^(n-1)*exp(-θ*f*modified_V)*√((θ*f*modified_V)/π) )/(θ*f)
    end
end

prefactor(θ,α,V,f,correction) = Real(sech(α))*exp(θ*(g(α)*θ^2*correction - V))/√f

function Z_integrand(θ,α,V,f,correction,modified_V)
    #Calculates the integrand of the partition function
    prefactor(θ,α,V,f,correction)*I(θ,f,modified_V,0)
end

function Z_integrand(q::Real,θ::Real,D::Real,Ω::Function,V::Function,correction::Function,par)
    #Calculates the integrand of the partition function
    α = θ*Ω(q,par)/2
    a_V = V(q,par)
    Z_integrand(θ,α,a_V,f(α),correction(q,par),D-a_V)
end

function H_integrand(θ,α,V,modified_V,correction,f)
    #Calculates the integrand of the energy
    prefactor(θ,α,V,f,correction)*(I(θ,f,modified_V,1)+I(θ,f,modified_V,0)*V)
end

function H_integrand(q::Real,θ::Real,D::Real,Ω::Function,V::Function,correction::Function,par)
    #Calculates the integrand of the energy
    α = θ*Ω(q,par)/2
    a_V = V(q,par)
    H_integrand(θ,α,a_V,D-a_V,correction(q,par),f(α))
end

function H2_integrand(θ,α,V,modified_V,correction,f)
    #Calculates the integrand of the energy
    prefactor(θ,α,V,f,correction)*(I(θ,f,modified_V,2)+2*I(θ,f,modified_V,1)*V + I(θ,f,modified_V,0)*V^2)
end

function H2_integrand(q::Real,θ::Real,D::Real,Ω::Function,V::Function,correction::Function,par)
    #Calculates the integrand of the energy
    α = θ*Ω(q,par)/2
    a_V = V(q,par)
    H2_integrand(θ,α,a_V,D-a_V,correction(q,par),f(α)) - real(Ω(q,par)/2)^2
end

function ode_step_energy(y,(θ,V,D,correction,Ω,par),q)
    dy1 = Z_integrand(q,θ,D,Ω,V,correction,par)
    dy2 = H_integrand(q,θ,D,Ω,V,correction,par)
    SA[dy1,dy2]
end

function ode_step_energy2(y,(θ,V,D,correction,Ω,par),q)
    dy1 = Z_integrand(q,θ,D,Ω,V,correction,par)
    dy2 = H2_integrand(q,θ,D,Ω,V,correction,par)
    SA[dy1,dy2]
end

function energy(θs::AbstractArray,par,V,D,correction,Ω,q₋,q₊)

    prob = ODEProblem(ode_step_energy,zeros(SVector{2,typeof(q₋)}),(q₋,q₊),(θs[1],V,D,correction,Ω,par))

    function prob_func(prob,i,repeat)
        remake(prob,p=(θs[i],V,D,correction,Ω,par))
    end

    output_func(sol,i) = (sol[2,end]/sol[1,end],false)

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)
    solve(ensemble_prob,trajectories=length(θs))
end

function energy(θs::AbstractArray,par,V,D,m,q₋,q₊)

    function correction(q,par)
        derivative(q->V(q,par),q)^2/(8m)
    end
    
    function Ω(q,par)
        √Complex(derivative(q->derivative(q->V(q,par),q),q)/m)
    end

    energy(θs,par,V,D,correction,Ω,q₋,q₊)
end

function energy2(θs::AbstractArray,par,V,D,correction,Ω,q₋,q₊)

    prob = ODEProblem(ode_step_energy2,zeros(SVector{2,typeof(q₋)}),(q₋,q₊),(θs[1],V,D,correction,Ω,par))

    function prob_func(prob,i,repeat)
        remake(prob,p=(θs[i],V,D,correction,Ω,par))
    end

    output_func(sol,i) = (sol[2,end]/sol[1,end],false)

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)
    solve(ensemble_prob,trajectories=length(θs))
end

function heat(θs::AbstractArray,par,V,D,m,q₋,q₊)

    function correction(q,par)
        derivative(q->V(q,par),q)^2/(8m)
    end
    
    function Ω(q,par)
        √Complex(derivative(q->derivative(q->V(q,par),q),q)/m)
    end

    (θs.^2).*(energy2(θs,par,V,D,correction,Ω,q₋,q₊) - energy(θs,par,V,D,correction,Ω,q₋,q₊).^2)
end