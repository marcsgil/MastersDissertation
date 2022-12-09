module SemiClassicalPöschlTeller

export sc_pt

__revise_mode__ = :eval

include("../../DoublePhaseSpace.jl")
include("../../ThermodynamicSystems.jl")

using LinearAlgebra,FastGaussQuadrature,ForwardDiff,StaticArrays
using .DoublePhaseSpace

V(q,yp,χ) = 0.5*(1/χ-χ)*(1+sech(2q)^2*(cos(yp)^2-2))/(1+sech(2q)*cos(yp))^2
fy(y,x,χ) = SA[-2χ*x[1], -ForwardDiff.derivative( q->V(q,y[1],χ),x[2] )]
fx(y,x,χ) = SA[ForwardDiff.derivative( yp->V(x[2],yp,χ),y[1] ), -0.5*χ*y[2]]
H(x,χ) = χ*x[1]^2+0.25*(1/χ-χ)*tanh(x[2])^2

function get_nws(N::Int)
    (ns_gl,ws_gl) = gausslegendre(N)
    (ns_gc,ws_gc) = gausschebyshev(N)
    Nodes(θ,χ) =  [ [0.5*√( (1/χ^2-1)*(1-Q^2) )*P,atanh(Q)] for P in ns_gl, Q in ns_gc ] 
    weights = [ w_gl*w_gc for w_gl in ws_gl, w_gc in ws_gc ]
    Nodes,weights
end

Nodes, weights = get_nws(60)

U(θ,χ) = DoublePhaseSpace.U(θ,χ,fy,fx,H,Nodes,weights)
C(θ,χ) = 0.0

sc_pt = ThermodynamicalSystem(U,C)
end