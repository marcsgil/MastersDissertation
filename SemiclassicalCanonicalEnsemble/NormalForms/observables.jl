include("wigner_symbols.jl")

function energy(θs::AbstractArray,H_operator::Polynomial,expectation_value::Function,J_max=200.)
    F = wigner_symbol_of_normal_form(H_operator)
    expectation_value( θs,F,F,J_max)
end

function heat_capacity(θs::AbstractArray,H_operator::Polynomial,expectation_value::Function,J_max=200.)
    F = wigner_symbol_of_normal_form(H_operator)
    F2 = wigner_symbol_of_normal_form(H_operator^2)
    Us = expectation_value( θs,F,F,J_max)
    Us2 = expectation_value( θs,F,F2,J_max)
    map( (U2,U,θ)-> θ^2*(U2-U^2), Us2,Us,θs )
end