module DMMorse

using ThreadedIterables

include("../../DMStdHamiltonian.jl")

V(q,χ) = (1-exp(-q))^2/(4χ)

U(θs::AbstractArray,χ) = energy(θs,χ,V,1/4χ,1/2χ,-log(2),20)

function U(θs::AbstractArray,χs::AbstractArray)
    result = tmap(χ->U(θs,χ),χs)
    reduce(hcat,result)
end

C(θs::AbstractArray,χ) = heat(θs,χ,V,1/4χ,1/2χ,-log(2),20)

function C(θs::AbstractArray,χs::AbstractArray)
    result = tmap(χ->C(θs,χ),χs)
    reduce(hcat,result)
end

end