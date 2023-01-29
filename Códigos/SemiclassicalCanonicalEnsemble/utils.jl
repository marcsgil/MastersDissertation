#Some random functions used throughout the code

function cartesian_product(x,y)
    result = fill((x[1],y[1]),length(x),length(y))
    for n in eachindex(y)
        for m in eachindex(x)
            result[m,n] = (x[m],y[n])
        end
    end
    result
end

function relative_error(x,y)
    min(abs(1-x/y),prevfloat(Inf))
end

get_last_half(v::AbstractArray) = @view v[length(v)÷2+1:end]

function twoD_nws((ns_1,ws_1),(ns_2,ws_2))
    [(ns_1[m],ns_2[n]) for m in eachindex(ns_1), n in eachindex(ns_2) ], 
    [ws_1[m]*ws_2[n] for m in eachindex(ws_1), n in eachindex(ws_2) ]
end