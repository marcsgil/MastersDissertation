using Plots,LaTeXStrings
theme(:wong2)
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)
##
V(q) = (1-exp(-q))^2
qs = LinRange(-1,4,600)
Vs = V.(qs)
plot(qs,Vs,
xlims=(first(qs),last(qs)),
ylims=(-.1,1.6),
label=false,
yticks = ([0,.5,1,1.5],[0,L"D/2",L"D",L"3D/2"]),
minorticks=false,
xlabel=L"a(r-r_e)",
ylabel=L"V(r)",
size=(700,400),
left_margin=2Plots.mm,
bottom_margin=2Plots.mm,color=:royalblue1
)
##
V(q,v₀) = v₀*(q^4-2q^2)
qs = LinRange(-2,2,600)
Vs = V.(qs,16)
plot(qs,Vs,
xlims=(first(qs),last(qs)),
ylims=(-18,10),
label=false,
yticks = ([-16,-8,0,8],[L"-V_0",L"-V_0/2",L"0",L"V_0/2"]),
minorticks=false,
xlabel=L"x/l",
ylabel=L"V(x)",
size=(700,400),
left_margin=2Plots.mm,
bottom_margin=2Plots.mm,color=:royalblue1
)