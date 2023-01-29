##
includet("Systems/Pöschl-Teller/ExactPöschlTeller.jl")
includet("Systems/Pöschl-Teller/DMPöschlTeller.jl")
includet("Systems/Pöschl-Teller/ClassicalPöschlTeller.jl")
includet("Systems/Pöschl-Teller/SemiClassicalPöschlTeller.jl")

using .ExactPöschlTeller,.SemiClassicalPöschlTeller,.DMPöschlTeller,.ClassicalPöschlTeller
##

##
using BenchmarkTools
using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5)
scalefontsizes(1.3)
##

##
function plot_relative_error!(error,grid)
    Threads.@threads for n in eachindex(grid)
        error[n] = min( prevfloat(Inf), 
        abs( 1-dm_pt.U(grid[n][1],grid[n][2])/ex_pt.U(grid[n][1],grid[n][2]) )  )
    end
    heatmap(τs,χs,error,clims=(0,0.1),dpi = 400)
end

function plot_functions(τs,functions,p;label=nothing)
    values = map(τ->functions[1](τ,p),τs)
    plots = plot(τs,values,xlims=(first(τs),last(τs)), ylims=(0,maximum(values)),
    dpi=400,label=label[1],xlabel=L"\beta \hbar \omega",
    ylabel=L"U/\hbar \omega",title=L"\chi=%$p")
    
    for n in 2:length(functions)
        plot!(τs,map(τ->functions[n](τ,p),τs),label=label[n])
    end

    plots
end
##

##
N = 100
τs = LinRange(0.4,5,N)
χs = LinRange(0.5,0.04,N)
χ = 0.005;
##

plot_functions(τs,[ex_pt.U,sc_pt.U,dm_pt.U,cl_pt.U],χ,label=["Exact" "SC" "DM" "Classical"])
savefig("Systems/Pöschl-Teller/Plots/χ=$χ.png")
plot_functions(τs,[ex_pt.U,dm_pt.U,cl_pt.U],χ,label=["Exact" "DM" "Classical"])

@benchmark sc_pt.U(1.0,0.025)

points = [ (τ,χ) for χ in χs, τ in τs ]
errors = zeros(size(points));

plot_relative_error!(errors, points)