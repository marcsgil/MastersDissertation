using FastGaussQuadrature,StaticArrays,LinearAlgebra,BenchmarkTools,ThreadedIterables,ComponentArrays,DifferentialEquations,SparseDiffTools,LoopVectorization
using Parameters: @unpack

using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5)
scalefontsizes(1.3)

Uex(θ,ω) = 0.5*ω*coth(0.5*θ*ω);

dy(y,x,ω) = -2*ω*x
dx(y,x,ω) = -0.5*ω*y
H(x,ω) = 0.5*ω*(x⋅x)

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

u0(x0) = ComponentArray( y=SA[0.0,0.0],x=x0,My= zeros(2,2),Mx=[1.0 0.0;0.0 1.0],Δ=0.0 )

function Z_integrand(u)
    exp( u.Δ + 0.5*log( abs(det(u.Mx)) )  )
end

coord_transformation(X,ω,θ) = √(2/(θ*ω))*X

function build_nws((xs1,ws1),(xs2,ws2))
    [ SA[xs1[m],xs2[n]] for m in eachindex(xs1), n in eachindex(xs2) ], [ws1[m]*ws2[n] for m in eachindex(ws1), n in eachindex(ws2) ]
end

f1(x) = x[1]
f2(x) = x[1]*x[2]

function Usc(θ,ω,nodes,weights)
    function prob_func(prob,i,repeat)
        remake(prob,u0=u0(coord_transformation(nodes[i],ω,θ)))
    end
    
    output_func = (sol,i) -> ([weights[i]*Z_integrand(last(sol.u)),H(last(sol.u).x,ω)],false)

    prob = ODEProblem( F!, u0(coord_transformation(nodes[1],ω,θ)), (0.0,0.5*θ), ω)
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)

    sols = solve(ensemble_prob,save_everystep=false,trajectories=length(nodes),reltol=1e-1,abstol=1e-2,verbose=false)

    sum(map(f2,sols))/(sum(map(f1,sols)))
end

##
nodes,weights = build_nws(gausshermite(100),gausshermite(100))

θ=5
ω=1

100*(Usc(θ,ω,nodes,weights)/Uex(θ,ω)-1)
##
Uex(θ,ω)
##
@benchmark Usc(θ,ω,nodes,weights)
##
θs = LinRange(0,5,50)
Ussc = tmap(θ->Usc(θ,ω,nodes,weights), θs)
Usex = tmap(θ->Uex(θ,ω),θs)
plot(θs,[Usex,Ussc])
##