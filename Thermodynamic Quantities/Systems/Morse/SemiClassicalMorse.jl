using ComponentArrays,DifferentialEquations,SparseDiffTools,LinearAlgebra,StaticArrays,ThreadedIterables,FastGaussQuadrature,BenchmarkTools,ColorSchemes,JLD
using Parameters: @unpack

using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)

includet("ExactMorse.jl")
includet("DampenedMetaplecticMorse.jl")
includet("ClassicalMorse.jl")

using .ExactMorse,.DampenedMetaplecticMorse,.ClassicalMorse

dy(y,x,χ) = SA[-4χ*x[1], ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ]
dx(y,x,χ) = SA[ (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ), -χ*y[2]]
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)

function Z_integrand(node,u,θ,χ)
    exp( u.Δ - θ*( node[1]^2 + node[2]^2 -(node[1]*node[2])^2 )/(4χ) + 0.5*log( abs(det(u.Mx)) )  )
end

coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]

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

function prob_func(prob,i,repeat,nodes,par)
    remake(prob,u0=u0(coord_transformation(nodes[i],par)))
end

function Usc(θs::AbstractArray,par,(nodes,weights),cb)

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

function build_nws((xs_gl,ws_gl),(xs_gc,ws_gc))
    [ (xs_gl[m],xs_gc[n]) for m in eachindex(xs_gl), n in eachindex(xs_gc) ], [ws_gl[m]*ws_gc[n] for m in eachindex(ws_gl), n in eachindex(ws_gc) ]
end
##
N=400
θs = LinRange(0.001,7,400)
χ= 0.08
nodes,weights = build_nws(gausslegendre(N), gausschebyshev(N,3))

Us = zeros(length(θs),4)
Us[:,1] = tmap(θ->ex_morse.U(θ,χ),θs)
Us[:,2] = Usc(θs,χ,(nodes,weights),cb)
Us[:,3] = tmap(θ->dm_morse.U(θ,χ),θs)
Us[:,4] = tmap(θ->cl_morse.U(θ,χ),θs)
##
default()
default(fontfamily="Computer Modern",
        linewidth=2.5, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.5)
p3 = plot(θs,Us,
xlims=(first(θs),last(θs)), ylims=(0,2.1),
xlabel = L"\theta", ylabel= L"U",
annotations = ((.5,.9),L"\chi=%$χ"), 
#labels = reshape( ["Exact","Semiclassical", "Dampened Metaplectic", "Classical"], 1, 4 ),
palette = [:black,:red,:blue,:orange],
line = reshape( [:solid, :dash, :dashdot, :dot],1,4 ),
left_margin=5Plots.mm,bottom_margin=5Plots.mm,label=false
)
##
plot(p1,p2,p3,p4,size=(1200,800),layout=(2,2))
##
@btime Usc(θs,χ,(nodes,weights),cb)
##
#@benchmark dm_morse.U(last(θs),last(χs))
relative_error(x,y) = min(prevfloat(Inf),abs(1-x/y))
χs = LinRange(0.005,0.12,400)

sc_results = load("Systems/Morse/Resultados/Erro Relativo Energia/error_sc.jld")
Ussc2 = sc_results["errors"]
θs = sc_results["θs"]
χs = sc_results["χs"]

Usex2 = [ex_morse.U(θ,χ) for χ in χs,θ in θs]
Uscl2 = [dm_morse.U(θ,χ) for χ in χs,θ in θs]
##
default()
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(2)
errors=tmap(relative_error,Uscl2,Usex2)
heatmap(θs,χs,Ussc2*100,clims=(0,8),xlabel=L"\theta",ylabel=L"\chi",yticks=[.02,.04,.06,.08,.1,.12],xticks=[0,1,2,3,4,5],size=(600,500))
##
savefig("error_cl.png")
##
χ = 0.08
prob = ODEProblem( F!, u0(SA[-0.06
-0.67]), (0.0,2), χ)
sol = solve( prob,Vern6(),reltol=1e-6)
##
function S(sol,χ,θs)
    map( θ-> sol(.5*θ).Δ-.5*θ*H(first(sol.u).x,χ),θs )
end

function M(sol,θs)
    map( θ-> det(sol(.5*θ).Mx),θs )
end
θs = LinRange(2.0,2.125,400)
Ss = S(sol,0.08,θs)
Ms = M(sol,θs)
##
plot(θs,Ss,label=false,right_margin=36Plots.mm,left_margin=3Plots.mm,xlims=(θs[begin],θs[end]),xlabel=L"\theta",ylabel=L"S_\theta^E",size=(650,400))
plot!(twinx(),θs,Ms,color=:red,xticks=:none,label=false,legend=:bottomright,xlims=(θs[begin],θs[end]),ylabel=L"\det \frac{\partial \mathbf{x}}{\partial\mathbf{X}}",line=:dash)
##
N=600
χ = .08
qs = LinRange(-log(2),2,N)
ps = LinRange(-1/(2χ),1/(2χ),N)
xs = [ SA[p,q] for p in ps, q in qs ]

function get_quantities(θ,χ,xs)

    function Z_integrand(node,u,θ,χ)
        exp( u.Δ - θ*H(node,χ) + 0.5*log( abs(det(u.Mx)) )  )
    end

    function prob_func(prob,i,repeat)
        remake(prob,u0=u0(xs[i]))
    end

    function output_func(sol,i)
		if H(xs[i],χ)>1/(4χ)
			([NaN,NaN],false)
		else
        	([Z_integrand(xs[i],last(sol.u),θ,χ),det(last(sol.u).Mx)],false)
		end
    end

    prob = ODEProblem( F!, u0(xs[1]), (0.0,0.5*θ), χ)
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func)

	sols = solve( ensemble_prob,Vern6(),trajectories=length(xs),verbose=false,reltol=1e-6,callback=cb )


	reshape(map(x->x[1],sols),(N,N)),reshape(map(x->x[2],sols),(N,N))
end
##
θ=3
Zs,dets = get_quantities(θ,χ,xs)

h1 = heatmap(qs,ps,Zs,clims=(0,1),xlabel=L"q",ylabel=L"p",title=L"\theta=3")
h2 = heatmap(qs,ps,dets,color=cgrad([:red,:white,:blue]),clims=(-2,2),xlabel=L"q",ylabel=L"p")
plot(h1,h2,size=(1200,480),plot_title = L"\theta=%$θ",top_margin=10Plots.mm,bottom_margin=6Plots.mm,left_margin=6Plots.mm)
savefig("det_part_3.png")