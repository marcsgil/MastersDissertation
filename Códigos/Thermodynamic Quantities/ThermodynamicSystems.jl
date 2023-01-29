using ForwardDiff:derivative

struct ThermodynamicalSystem
    U::Function #Energy
    C::Function #Specific Heat
end

function ThermodynamicalSystem(U)
    ThermodynamicalSystem(U, (τ,p)->-τ^2*derivative( τ-> U(τ, p), τ ))
end

function thermo_quantities_from_partition(Z)
    U(τ,p) = -derivative( τ->log(Z(τ,p)), τ )
    C(τ,p)= -τ^2*derivative( τ-> U(τ, p), τ )
    U,C
end

function partition_from_spectrum( spectrum::Function, τ, p, N::Int64, tol::Float64 )
    sum = exp( -τ*spectrum(0,p) )
    for n in 1:N
        a = exp(-τ*spectrum(n,p))
        sum +=a
        if a/sum < tol
            break
        else
        end
    end
    sum
end

function thermo_quantities_from_spectrum(spectrum::Function, N, tol)
    thermo_quantities_from_partition((τ,p)->partition_from_spectrum(spectrum,τ,p,N,tol))
end