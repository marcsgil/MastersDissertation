using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
#scalefontsizes(1.3)
scalefontsizes(1.5)
##
using Revise
includet("EXMorse.jl")
includet("SCMorse.jl")
includet("DMMorse.jl")
includet("CLMorse.jl")
includet("../../utils.jl")
##
θs = LinRange(0.001,7,300)
χ = 5e-3
##
Usex = EXMorse.U.(θs,χ)
Ussc = SCMorse.U(θs,χ,400)
Usdm = DMMorse.U(θs,χ)
Uscl = CLMorse.U.(θs,χ)
##
plot(θs,[Usex,Ussc,Usdm,Uscl],
xlims=(first(θs),last(θs)), ylims=(0,4))
#=plot(θs,[Usex,Ussc],
xlims=(first(θs),last(θs)), ylims=(0,4))=#
##
@benchmark SCMorse.U($θs,$χ,150)
χs = LinRange(5e-3,.12,300)
Ussc = SCMorse.U(θs,χs,150)
@benchmark SCMorse.U($θs,$χs,130)
##
pars = cartesian_product(θs,χs)
Usex = map(par->EXMorse.U(par[1],par[2]),pars)
error = map(relative_error,Ussc,Usex)*100
##
heatmap(θs,χs,error',clims=(0,8),yticks=[.02,.04,.06,.08,.1,.12],xlabel=L"\theta",ylabel=L"\chi")
plot!(4χs,χs,label=false,color=:green,xlims)
##
θs = LinRange(.001,3,300)
χ = .005
Csex = EXMorse.C.(θs,χ)
Cssc = SCMorse.C(θs,χ,80)
Csdm = DMMorse.C(θs,χ)
Cscl = CLMorse.C.(θs,χ)
plot(θs,[Csex,Cssc,Csdm,Cscl],ylims=(0,1.5))
##
@benchmark SCMorse.C($θs,$χ,80)
##
θs = LinRange(.001,7,100)
function make_plot(χ)
        Us = zeros(length(θs),4)
        Us[:,1] = EXMorse.U.(θs,χ)
        Us[:,2] = SCMorse.U(θs,χ,300)
        Us[:,3] = DMMorse.U(θs,χ)
        Us[:,4] = CLMorse.U.(θs,χ)
    
        ylims=(0,min(5*Us[end,1],Us[begin,1]+.1))
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
χs = [.01,.04,.08,.12]
plot(ntuple(n->make_plot(χs[n]),4)...,layout=(2,2),size=(1200,800))
##
θs = LinRange(.001,5,300)
function make_plot(χ)
        Cs = zeros(length(θs),4)
        Cs[:,1] = EXMorse.C.(θs,χ)
        Cs[:,2] = SCMorse.C(θs,χ,200)
        Cs[:,3] = DMMorse.C(θs,χ)
        Cs[:,4] = CLMorse.C.(θs,χ)
    
        ylims=(0,maximum(Cs[:,4])*1.1)
        plot(θs,Cs,
        xlims=(first(θs),last(θs)), ylims=ylims,
        xlabel = L"\theta", ylabel= L"c",
        annotations = ((.5,.9),L"\chi=%$χ"), 
        #labels = reshape( ["Exact","Semiclassical", "Dampened Metaplectic", "Classical"], 1, 4 ),
        palette = [:black,:red,:blue,:orange],
        line = reshape( [:solid, :dash, :dashdot, :dot],1,4 ),
        left_margin=5Plots.mm,bottom_margin=5Plots.mm,label=false
        )
end
χs = [.01,.04,.08,.12]
plot(ntuple(n->make_plot(χs[n]),4)...,layout=(2,2),size=(1200,800))
#savefig("Systems/Morse/Resultados/calores_morse.png")
##
N = 300
θs = LinRange(10^-3,2,N)
χs = LinRange(.005,.12,N)
pars = cartesian_product(θs,χs)
Csex = map(par->EXMorse.C(par[1],par[2]),pars)
##
Cssc = SCMorse.C(θs,χs,80)
error = map(relative_error,Cssc,Csex)*100
h_usc = heatmap(θs,χs,error',
xlims=(0,last(θs)), ylims=(first(χs),last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.02, .04 ,.06, .08 ,.1,.12])
##
Csdm = DMMorse.C(θs,χs)
error = map(relative_error,Csdm,Csex)*100
h_usc = heatmap(θs,χs,error',
xlims=(0,last(θs)), ylims=(first(χs),last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.02, .04 ,.06, .08 ,.1,.12])
##
Cscl = map(par->CLMorse.C(par[1],par[2]),pars)
error = map(relative_error,Cscl,Csex)*100
h_usc = heatmap(θs,χs,error',
xlims=(0,last(θs)), ylims=(first(χs),last(χs)),clims=(0,8),
xlabel = L"\theta", ylabel= L"χ",size=(600,500), yticks = [.02, .04 ,.06, .08 ,.1,.12])