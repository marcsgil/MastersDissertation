module SCKerr

include("../../NormalForms/SCNormalForms.jl")
include("../../NormalForms/observables.jl")

U(θs::AbstractArray,χ) = energy(θs,Polynomial([0,1,χ]),expectation_value,2000.)
C(θs::AbstractArray,χ) = heat_capacity(θs,Polynomial([0,1,χ]),expectation_value,2000.)

function U(θs::AbstractArray,χs::AbstractArray)
    reduce(hcat,map(χ->U(θs,χ),χs))
end

function C(θs::AbstractArray,χs::AbstractArray)
    reduce(hcat,map(χ->C(θs,χ),χs))
end

end