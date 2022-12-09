using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
#scalefontsizes(1.3)
scalefontsizes(2)
##
using Revise
includet("EXKerr.jl")
includet("SCKerr.jl")
includet("DMKerr.jl")
includet("CLKerr.jl")
includet("../../utils.jl")
##
θs = LinRange(.1,12,100)
χ = .002
Usex = EXKerr.U.(θs,χ)
Ussc = SCKerr.U(θs,χ)
Usdm = DMKerr.U(θs,χ)
Uscl = CLKerr.U.(θs,χ)
plot(θs,[Usex,Ussc,Usdm,Uscl],ylims=(.4,5),xlims=(first(θs),last(θs)))
##
θs = LinRange(.1,12,100)
function make_plot(χ)
        Us = zeros(length(θs),4)
        Us[:,1] = EXKerr.U.(θs,χ)
        Us[:,2] = SCKerr.U(θs,χ)
        Us[:,3] = DMKerr.U(θs,χ)
        Us[:,4] = CLKerr.U.(θs,χ)
    
        ylims=(0,2.5*Us[end,1])
        plot(θs,Us,
        xlims=(first(θs),last(θs)), ylims=ylims,
        xlabel = L"\theta", ylabel= L"U",
        annotations = ((.5,.9),L"\chi=%$χ"), 
        #labels = reshape( ["Exact","Semiclassical", "Dampened Metaplectic", "Classical"], 1, 4 ),
        palette = [:black,:red,:blue,:orange],
        line = reshape( [:solid, :dash, :dashdot, :dot],1,4 ),
        left_margin=5Plots.mm,bottom_margin=5Plots.mm,label=false
        )
end
χs = [.01,.1,.5,1.]
plot(ntuple(n->make_plot(χs[n]),4)...,layout=(2,2),size=(1200,800))
savefig("Systems/Kerr/Resultados/energias_kerr.png")
##
N = 300
θs = LinRange(10^-3,2,N)
χs = LinRange(10^-3,1,N)
pars = cartesian_product(θs,χs)
Usex = map(par->EXKerr.U(par[1],par[2]),pars)
##
Ussc = SCKerr.U(θs,χs)
error_sc = map(relative_error,Ussc,Usex)*100
h_usc = heatmap(θs,χs,error_sc',
xlims=(0,last(θs)), ylims=(first(χs),last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
Usdm = DMKerr.U(θs,χs)
error_dm = map(relative_error,Usdm,Usex)*100
h_udm = heatmap(θs,χs,error_dm',
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
Uscl = map(par->CLKerr.U(par[1],par[2]),pars)
error_cl = map(relative_error,Uscl,Usex)*100
h_ucl = heatmap(θs,χs,error_cl',
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
θs = LinRange(.1,2,100)
χ = .1
Csex = EXKerr.C.(θs,χ)
Cssc = SCKerr.C(θs,χ)
Csdm = DMKerr.C(θs,χ)
Cscl = CLKerr.C.(θs,χ)
plot(θs,[Csex,Cssc,Csdm,Cscl],ylims=(0,1))
##
θs = LinRange(.1,10,100)
function make_plot(χ)
        Cs = zeros(length(θs),4)
        Cs[:,1] = EXKerr.C.(θs,χ)
        Cs[:,2] = SCKerr.C(θs,χ)
        Cs[:,3] = DMKerr.C(θs,χ)
        Cs[:,4] = CLKerr.C.(θs,χ)
    
        ylims=(-.05,1.1*Cs[end,4])
        plot(θs,Cs,
        xlims=(first(θs),last(θs)), ylims=ylims,
        xlabel = L"\theta", ylabel= L"c",
        annotations = ((.8,.7),L"\chi=%$χ"), 
        #labels = reshape( ["Exact","Semiclassical", "Dampened Metaplectic", "Classical"], 1, 4 ),
        palette = [:black,:red,:blue,:orange],
        line = reshape( [:solid, :dash, :dashdot, :dot],1,4 ),
        left_margin=5Plots.mm,bottom_margin=5Plots.mm,label=false
        )
end
χs = [.01,.1,.5,1.]
plot(ntuple(n->make_plot(χs[n]),4)...,layout=(2,2),size=(1200,800))
##
N = 300
θs = LinRange(10^-3,2,N)
χs = LinRange(10^-3,1,N)
pars = cartesian_product(θs,χs)
Csex = map(par->EXKerr.C(par[1],par[2]),pars)
##
Cssc = SCKerr.C(θs,χs)
error_sc = map(relative_error,Cssc,Csex)*100
h_csc = heatmap(θs,χs,error_sc',
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
Csdm = DMKerr.C(θs,χs)
error_dm = map(relative_error,Csdm,Csex)*100
h_cdm = heatmap(θs,χs,error_dm',
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##
Cscl = map(par->CLKerr.C(par[1],par[2]),pars)
error_cl = map(relative_error,Cscl,Csex)*100
h_ccl = heatmap(θs,χs,error_cl',
xlims=(0,last(θs)), ylims=(0,last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.2, .4 ,.6, .8 ,1])
##