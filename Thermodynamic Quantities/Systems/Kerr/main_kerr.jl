using Revise
includet("ExactKerr.jl")
includet("SemiClassicalKerr.jl")
includet("DampenedMetaplecticKerr.jl")
includet("ClassicalKerr.jl")

using BenchmarkTools,.ExactKerr,.SemiClassicalKerr,.DampenedMetaplecticKerr,.ClassicalKerr,ThreadedIterables

using Plots,LaTeXStrings,ColorSchemes


relative_error(x,y) = min(prevfloat(Inf),abs(1-x/y))
relative_error(f,g,grid) = tmap( x -> relative_error(f(x...),g(x...)),grid)
##
#Parameters
N = 100
θs = LinRange(10^-3,2,N)

#errors = relative_error(sc_kerr.U,ex_kerr.U,points)
#heatmap(θs,χs,errors,color=:turbo,clims=(0,0.1),xlabel=L"\theta",ylabel=L"\chi")
##
@benchmark sc_kerr.U(1,1)
##
theme(:default)
default(fontfamily="Computer Modern",
        linewidth=2.5, framestyle=:box, label=nothing, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)
##
function make_plot(χ)
    Us = zeros(length(θs),4)
    Us[:,1] = tmap(θ->ex_kerr.U(θ,χ),θs)
    Us[:,2] = tmap(θ->sc_kerr.U(θ,χ),θs)
    Us[:,3] = tmap(θ->dm_kerr.U(θ,χ),θs)
    Us[:,4] = tmap(θ->cl_kerr.U(θ,χ),θs)

    #ylims=(0,2.5*Us[end,1])
    plot(θs,Us,
    xlims=(first(θs),last(θs)), ylims=(0,1),
    xlabel = L"\theta", ylabel= L"U",
    annotations = ((.5,.9),L"\chi=%$χ"), 
    #labels = reshape( ["Exact","Semiclassical", "Dampened Metaplectic", "Classical"], 1, 4 ),
    palette = [:black,:red,:blue,:orange],
    line = reshape( [:solid, :dash, :dashdot, :dot],1,4 ),
    left_margin=5Plots.mm,bottom_margin=5Plots.mm
    )
end

χs = [.01,.1,.5,1.]
plot(ntuple(n->make_plot(χs[n]),4)...,layout=(2,2),size=(1200,800))
##
theme(:default)
default(fontfamily="Computer Modern",
        linewidth=2.5, framestyle=:box, label=nothing, grid=true, minorticks=5,dpi=400)
scalefontsizes(2)
N = 600
θs = LinRange(10^-3,2,N)
χs = LinRange(10^-3,1,N)
points = [ (θ,χ) for χ in χs, θ in θs ]
error_sc = relative_error(sc_kerr.U,ex_kerr.U,points)*100
error_dm = relative_error(dm_kerr.U,ex_kerr.U,points)*100
error_cl = relative_error(cl_kerr.U,ex_kerr.U,points)*100
##
h_sc = heatmap(θs,χs,error_sc,
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
#
h_dm = heatmap(θs,χs,error_dm,
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
#
h_cl = heatmap(θs,χs,error_cl,
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
theme(:default)
default(fontfamily="Computer Modern",
        linewidth=2.5, framestyle=:box, label=nothing, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.5)
N = 600
θs = LinRange(10^-3,12,N)
χs = LinRange(10^-3,1,N)
points = [ (θ,χ) for χ in χs, θ in θs ]
error_sc = relative_error(sc_kerr.U,ex_kerr.U,points)*100
error_dm = relative_error(dm_kerr.U,ex_kerr.U,points)*100
error_cl = relative_error(cl_kerr.U,ex_kerr.U,points)*100
##
h_sc = heatmap(θs,χs,error_sc,
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ", yticks = [.2, .4 ,.6, .8 ,1])