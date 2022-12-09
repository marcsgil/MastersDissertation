using ForwardDiff:derivative

function partition_from_spectrum( spectrum::Function, θ, p, N::Int64, tol::Float64 )
    sum = exp( -θ*spectrum(0,p) )
    for n in 1:N
        a = exp(-θ*spectrum(n,p))
        sum +=a
        if a/sum < tol
            break
        else
        end
    end
    sum
end

function energy_from_partition(Z::Function)
    (θ,par) -> -derivative( θ->log(Z(θ,par)), θ )
end

function heat_from_partition(Z::Function)
    (θ,par) -> θ^2*derivative(θ->derivative( θ->log(Z(θ,par)), θ ),θ)
end