using Plots,LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)
##
fr(x) = ifelse(iszero(x),1,tanh(x)/x)
fi(x) = ifelse(iszero(x),1,real(fr(im*x)))
gr(x) = real((1-fr(x))/x^2)
gi(x) = real(-(1-fi(x))/x^2)
##
xs = LinRange(-10,10,400)
plot(xs,fr.(xs),xlims=(first(xs),last(xs)),label=false,xlabel=L"\alpha",yticks=([0,1/3,2/3,1],["0","1/3","2/3","1"]),ylims=(0,1.05),size=(500,400))
plot!(xs,gr.(xs),xlims=(first(xs),last(xs)),label=false,xlabel=L"\alpha",yticks=([0,1/3,2/3,1],["0","1/3","2/3","1"]),ylims=(0,1.05),size=(500,400),color=:red,line=:dash)
savefig("fg.png")
##
xs = LinRange(-π/2+.01,π/2-.01,400)
plot(xs,fi.(xs),xlims=(first(xs),last(xs)),label=false,xlabel=L"\alpha",ylims=(0.2,4),size=(500,400))
plot!(xs,gi.(xs),xlims=(first(xs),last(xs)),label=false,xlabel=L"\alpha",size=(500,400),color=:red,line=:dash)
savefig("fgi.png")