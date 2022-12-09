module DoublePhaseSpace

export Usc,get_solution,u0,F!

using ComponentArrays,DifferentialEquations,SparseDiffTools,LinearAlgebra,StaticArrays,ThreadedIterables
using Parameters: @unpack


function F!(du,u,par,t)
    @unpack x,y,Mx,My,Δ = u

    du.y = dy(y,x,par)
    du.x = dx(y,x,par)

    for n in 1:size(My)[2]
        du.My[:,n] = auto_jacvec(y->dy(y,x,par), y, My[:,n]) + auto_jacvec(x->dy(y,x,par), x, Mx[:,n])
        du.Mx[:,n] = auto_jacvec(y->dx(y,x,par), y, My[:,n]) + auto_jacvec(x->dx(y,x,par), x, Mx[:,n])
    end

    du.Δ = y⋅dx(y,x,par)
end

u0(x0) = ComponentArray( y=SA[0.0,0.0],x=x0,My= zeros(2,2),Mx=[1.0 0.0;0.0 1.0],Δ=0.0)

condition(u,t,integrator) = det(u.Mx)
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

function Usc(θs::AbstractArray,par,(nodes,weights);alg=Tsit5(),reltol=1e-3,abstol=1e-6,callback=cb)

    prob = ODEProblem( F!, u0(coord_transformation(nodes[1],par)), (0.0,0.5*last(θs)), par)

    ensemble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat)->prob_func(prob,i,repeat,nodes,par))

    sols = solve( ensemble_prob,BS3(),trajectories=length(nodes),verbose=false,callback=cb,reltol=1e-1,abstol=1e-2 )

    result = similar(θs)
    Zs = similar(weights)
    Hs = similar(weights)

    for (i,θ) in enumerate(θs)
        map!( (node,weight,sol) -> weight*Z_integrand(node,sol(0.5*θ),θ,par), Zs,nodes,weights,sols )
        map!( (Z,sol) -> Z*H(sol(0.5*θ).x,par), Hs, Zs, sols )
        result[i] = sum(Hs)/sum(Zs)
    end
    
    result
end

function Usc(θs::AbstractArray,χs::AbstractArray,(nodes,weights),cb)

    result = zeros(length(χs),length(θs))

    Threads.@threads for k in eachindex(χs)
        result[k,:] = Usc(θs,χs[k],(nodes,weights),cb)
    end

    result
end

end



