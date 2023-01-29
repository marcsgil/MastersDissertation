using DiffEqSensitivity, DifferentialEquations, ForwardDiff

function lotka_volterra!(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end
p = [1.5,1.0,3.0,1.0]; u0 = [1.0;1.0]
prob = ODEProblem(lotka_volterra!,u0,(0.0,10.0),p)
sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6)


function f(x)
    _prob = remake(prob,u0=x[1:2],p=x[3:end])
    solve(_prob,Tsit5(),reltol=1e-6,abstol=1e-6,saveat=1)[1,:]
end
x = [u0;p]
f(x)
dx = ForwardDiff.jacobian(f,x)