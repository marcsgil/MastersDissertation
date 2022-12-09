module Utils

export  find_cutoff, evaluate_integral,evaluate_complex_integral,pmap,evaluate_integral_v2

using LinearAlgebra:⋅
using LoopVectorization

function evaluate_integral(integrand::Function, p, (nodes,weights)::Tuple{Vector{Float64}, Vector{Float64}})
    vmap( x-> integrand(x,p), nodes )⋅weights
end

function evaluate_integral_v2(integrand::Function, p, (nodes,weights)::Tuple{Vector{Float64}, Vector{Float64}})
    integrands = fill(integrand(first(nodes),p),size(nodes))
    vmapntt!( x-> integrand(x,p), integrands, nodes )
    integrands⋅weights
end

function pmap(f,collection)
    result = fill(f(collection[1]),size(collection))
    Threads.@threads for n in 2:length(collection)
        result[n] = f(collection[n])
    end
    result
end

function evaluate_complex_integral(integrand::Function, p, (nodes,weights)::Tuple{Vector{Float64}, Vector{Float64}})
    integrands = zeros(size(nodes))

    Threads.@threads for n in eachindex(integrands)
        integrands[n] = integrand(nodes[n],p)
    end
    
    integrands⋅weights
end

function find_cutoff(f,(tol,a,b,ϵ))
    #=Finds an interval of size I such that |I| ≤ ϵ 
    such that there is x ∈ I with f(x) < tol. 
    A bisection method is used, starting from the interval (a,b)=#

    for n in 1:ceil( Int, log(2,(b-a)/ϵ) )
        if abs( f(0.5*(a+b)) ) ≥ tol
            a=0.5*(a+b)
        else
            b=0.5*(a+b)
        end
    end
    return 0.5*(a+b)
end




end