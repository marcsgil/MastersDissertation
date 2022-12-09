using ThreadedIterables

function Z_integrand(u,energy_term)
    exp( u.Δ - energy_term )*√abs(det(u.jac_x))
end

#=function abstract_Usc(θs,H,dy,dx,(nodes,weights),par=nothing)

    sols = get_solutions(last(θs)/2,dy,dx,nodes,par)

    result = similar(θs)
    Zs = similar(weights)
    Hs = similar(weights)

    for (i,θ) in enumerate(θs)
        tmap!( (weight,sol) -> weight*Z_integrand(sol(0.5*θ),θ*H(sol(0).x,par)), Zs,weights,sols )
        tmap!( (Z,sol) -> Z*H(sol(0.5*θ).x,par), Hs, Zs, sols )
        result[i] = sum(Hs)/sum(Zs)
    end
    
    result
end

function abstract_Csc(θs,H,H2,dy,dx,(nodes,weights),par=nothing)

    sols = get_solutions(last(θs)/2,dy,dx,nodes,par)

    result = similar(θs)
    result2 = similar(θs)
    Zs = similar(weights)
    Hs = similar(weights)
    H2s = similar(weights)

    for (i,θ) in enumerate(θs)
        tmap!( (weight,sol) -> weight*Z_integrand(sol(0.5*θ),θ*H(sol(0).x,par)), Zs,weights,sols )
        tmap!( (Z,sol) -> Z*H(sol(0.5*θ).x,par), Hs, Zs, sols )
        tmap!( (Z,sol) -> Z*H2(sol(0.5*θ).x,par), H2s, Zs, sols )
        result[i] = sum(Hs)/sum(Zs)
        result2[i] = sum(H2s)/sum(Zs)
    end
    
    (θs.^2).*(result2.-result.^2)
end=#

function energy_output(sol,i,(nodes,weights),θs,par,H)
    output = fill( ntuple(i->zero(eltype(θs)),2), length(θs) )
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[n] = (weights[i]*Z_integrand(u,θs[n]*E),H(u.x,par))
    end
    output,false
end

function energy_reduction(sols,θs)
    f(collection) = sum( x->x[1]*x[2],collection)/sum( x->x[1],collection)
    tmap(n->f(view(sols,n,:)),eachindex(θs))
end

function heat_output(sol,i,(nodes,weights),θs,par,H,H2)
    output = fill( ntuple(i->zero(eltype(θs)),3), length(θs) )
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[n] = (weights[i]*Z_integrand(u,θs[n]*E),H(u.x,par),H2(u.x,par))
    end
    output,false
end

function heat_reduction(sols,θs)
    function f(collection)
        inverse_Z = 1/sum( x->x[1],collection)
        inverse_Z*sum( x->x[1]*x[3],collection)-(inverse_Z*sum( x->x[1]*x[2],collection))^2
    end
    tmap(n->θs[n]^2*f(view(sols,n,:)),eachindex(θs))
end