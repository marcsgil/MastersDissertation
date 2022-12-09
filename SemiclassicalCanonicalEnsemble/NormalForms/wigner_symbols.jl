using Polynomials

function get_converted_monomials(N)
    # monomials[n] =  ∑ cⱼ (p²+q²)ʲ, the Wigner representation of (̂p²+̂q²)ⁿ
    p = Polynomial([0,1])
    monomials = [p for n in 1:N+1]
    monomials[1] = one(p)
    for n in 3:N+1
        monomials[n] = p*monomials[n-1]-derivative(monomials[n-1])-p*derivative(derivative(monomials[n-1]))
    end
    monomials
end

function normalize_and_change_basis(p::Polynomial)
    #Given a basis 1,…,xⁿ and a polynomial p(x) = ∑ cₙ xⁿ, returns returns p/(2^degrre(p)) with respect to the basis 1,…,uⁿ, where u = x/2

    #This is useful because it is straightforward to obtain the Wigner representation of (̂p²+̂q²)ⁿ, which is done by get_converted_monomials.
    #On the other hand, to perform other calculations, it is easier to express normal forms as powers of [(̂p²+̂q²)/2]ⁿ or [(p²+q²)/2]ⁿ

    ds = coeffs(p)/1
    n = degree(p)
    for j ∈ eachindex(ds)
        ds[j]/=2^(n-j+1)
    end
    Polynomial(ds)
end

function wigner_symbol_of_normal_form(quantum_nf::Polynomial;assume_divided_by_two_basis=true)
    #Returns the Wigner symbol of a quantum normal form
    #If assume_divided_by_two_basis=true, we have quantum_nf = ∑ cₙ ̂uⁿ and result = ∑ dₙ uⁿ, where ̂u = (̂p²+̂q²)/2 and u = (p²+q²)/2
    #If assume_divided_by_two_basis=false, we have quantum_nf = ∑ cₙ ̂uⁿ and result = ∑ dₙ uⁿ, where ̂u = p²+̂q² and u = p²+q²

    converted_monomials = get_converted_monomials(degree(quantum_nf))

    if assume_divided_by_two_basis
        result = zero(quantum_nf)/1
        for (n,c) in enumerate(coeffs(quantum_nf))
            result += c*normalize_and_change_basis(converted_monomials[n])
        end
    else
        result = zero(quantum_nf)
        for (n,c) in enumerate(coeffs(quantum_nf))
            result += c*converted_monomials[n]
        end
    end

    result
end