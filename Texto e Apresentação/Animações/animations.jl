using DifferentialEquations, DiffEqPhysics, CairoMakie
##
H(p,q,ω) = p^2/2 + (1-exp(-q))^2

ts = LinRange(0,10π,360)
N = 19
qs = zeros(N)
ps = LinRange(-2,2,N)
prob_func(prob,i,repeat) = remake(prob,u0 =ArrayPartition((ps[i],qs[i])))
prob = HamiltonianProblem(H, ps[1], qs[1], (first(ts),last(ts)), 1.)
ensembleprob = EnsembleProblem(prob,prob_func=prob_func)
sols = solve(ensembleprob,trajectories = N)
##
fig = Figure(resolution=(800,800),fontsize=32)
ax = Axis(fig[1,1],xlabel="q",ylabel="p")
xlims!(ax,-1.1,3.1)
ylims!(ax,-2.1,2.1)

function update(t)
    ax.title = "t=$(round(t,digits=1))"
    empty!(ax)
    ps,qs = DifferentialEquations.EnsembleAnalysis.componentwise_vectors_timepoint(sols,t)
    scatter!(ax,qs,ps,color=:red,markersize=18)
    trail = LinRange(max(0,t-1),t,128)
    for u in sols
        lines!(ax,[ u(s)[2] for s in trail ],[ u(s)[1] for s in trail ],color = :blue)
    end
    fig
end
##
record(update,fig, "test.mp4",ts,fps=24)