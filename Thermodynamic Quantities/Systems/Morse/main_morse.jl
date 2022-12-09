using FastGaussQuadrature,StaticArrays,LinearAlgebra,ThreadedIterables,DifferentialEquations

using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)

includet("ExactMorse.jl")
includet("DampenedMetaplecticMorse.jl")
includet("ClassicalMorse.jl")

using .ExactMorse,.DampenedMetaplecticMorse,.ClassicalMorse

includet("../../DoublePhaseSpace.jl")
using .DoublePhaseSpace
using BenchmarkTools
using JLD


dy(y,x,χ) = SA[-4χ*x[1], ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ]
dx(y,x,χ) = SA[ (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ), -χ*y[2]]
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)

function Z_integrand(node,u,θ,χ)
    exp( u.Δ - θ*( node[1]^2 + node[2]^2 -(node[1]*node[2])^2 )/(4χ) + 0.5*log( abs(det(u.Mx)) )  )
end

coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]

function build_nws((xs_gl,ws_gl),(xs_gc,ws_gc))
    [ (xs_gl[m],xs_gc[n]) for m in eachindex(xs_gl), n in eachindex(xs_gc) ], [ws_gl[m]*ws_gc[n] for m in eachindex(ws_gl), n in eachindex(ws_gc) ]
end

relative_error(x,y) = min(prevfloat(Inf),abs(1-x/y))
relative_error(f,g,grid) = tmap( x -> relative_error(f(x),g(x)),grid)

function prob_func(prob,i,repeat,nodes)
    remake(prob,u0=u0(coord_transformation(nodes[i],prob.p)))
end

function output_func(sol,i,θ,χ,nodes,weights)
    ([weights[i]*Z_integrand(nodes[i],last(sol.u),θ,χ),H(last(sol.u).x,χ)],false)
end

function get_result(sols)
    #=N = Int(√length(sols))
    Zs = reshape(map(x->x[1],sols),(N,N))
    grad_x,grad_y = imgradients((@view Zs[N÷2:N÷2+1,:]), KernelFactors.ando3, "replicate")
    gradient = tmap((x,y)->x^2+y^2,grad_x,grad_y)

    j0=findlast(x-> x>20 || isnan(x) || isinf(x),@view gradient[1,:]) 
        
    for i in 1:N
        for j in 1:N
            Zs[i,j] = ifelse(j>ifelse(j0 == nothing,0,j0),Zs[i,j],0)
        end
    end

    sum( map((Z,y)->Z*y[2],vec(Zs),sols) )/sum(Zs)=#
    sum( map(x->x[1]*x[2], sols) )/sum(map(x->x[1],sols))
end
##
N = 30
nodes,weights = build_nws(gausslegendre(N),gausschebyshev(N))

θs = LinRange(0.001,2,100)
χ = 0.12

Usex = tmap(θ->ex_morse.U(θ,χ),θs)
Usdm = tmap(θ->dm_morse.U(θ,χ),θs)
Uscl = tmap(θ->cl_morse.U(θ,χ),θs)
Ussc = tmap(θ->Usc(θ,χ,(nodes,weights),dy,dx,prob_func,output_func,get_result),θs)

@benchmark tmap(θ->dm_morse.U(θ,χ),θs)


plot( θs,Usex,color=:black,
label="Exact",
xlabel=L"\theta",
ylabel=L"U",
xlims=(first(θs),last(θs)),
ylims=(0.5*last(Usex),min(first(Usex),4*last(Usex))),
annotations=((.5,.95),L"\chi=%$χ") )

plot!(θs,Ussc,label="Semi Classical",linestyle=:dot,color=:red,linewidth=3)
plot!(θs,Usdm,label="DM",linestyle=:dash,color=:blue)
plot!(θs,Uscl,label="Classical",linestyle=:dashdot,color=:green)
##
big_θs = LinRange(0.01,4,300)
big_χs = LinRange(0.01,0.1,300)
##
i_χ_min = 1
i_χ_max = 79
i_θ_min = 151
i_θ_max = 300
θs = big_θs[i_θ_min:i_θ_max]
χs = big_χs[i_χ_min:i_χ_max]
#θs = LinRange(big_θs[i_θ_min],big_θs[i_θ_max],20)
#χs = LinRange(big_χs[i_χ_min],big_χs[i_χ_max],20)
N = 100
nodes,weights = build_nws(gausslegendre(N),gausschebyshev(N))
relative_error(Usc(last(θs),first(χs),(nodes,weights),dy,dx,prob_func,output_func,get_result),ex_morse.U(last(θs),first(χs)))

@benchmark Usc(last(θs),last(χs),(nodes,weights),dy,dx,prob_func,output_func,get_result)
##
pars = [ (θ,χ) for χ in χs, θ in θs ]

Ussc = tmap( par -> Usc(par[1],par[2],(nodes,weights),dy,dx,prob_func,output_func,get_result),pars)
Usex = tmap( par -> ex_morse.U(par[1],par[2]), pars )


#errors = tmap( relative_error, Ussc,Usex )
errors = zeros(size(pars))
heatmap(θs,χs,errors,clims=(0,0.1),xlabel=L"\theta",ylabel=L"\chi")
##
save("Systems/Morse/relative_errors/quadrant3.jld","i_χ_min",i_χ_min,"i_χ_max",i_χ_max,"i_θ_min",i_θ_min,"i_θ_max",i_θ_max,"errors",errors)
errors1 = load("Systems/Morse/relative_errors/quadrant1.jld")["errors"]
errors2 = load("Systems/Morse/relative_errors/quadrant2.jld")["errors"]
errors3 = load("Systems/Morse/relative_errors/quadrant3.jld")["errors"]

heatmap(big_θs,big_χs,hcat(errors1,vcat(errors3,errors2)),clims=(0,0.1),xlabel=L"\theta",ylabel=L"\chi")
savefig("Systems/Morse/Novos/erro_relativo_sc2.png")