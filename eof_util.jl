
# function reshape_data(d::ncData)
#     M = length(d.lonvec) #size(d.data, 1)
#     N =  length(d.latvec)  #size(d.data, 2)
#     L = length(d.timevec)  # size(d.data, 3)
#     return reshape(d.data, (M * N, L))
# end

function reshape_data(data)
    M = size(data, 1)
    N = size(data, 2)
    L = size(data, 3)
    return reshape(data, (M * N, L))
end

function shape_data(data, M, N)
    return reshape(data, (M, N))
end

function get_eof(d::ncData, full=false)
    X = reshape_data(d.data)
    U, S, V = svd(X) 
    return full ? U : U,S,V
end


function projection(v, U)
    #for an orthogonal basis
    dim = size(U,2)
    proj = zeros(dim)
    for i in 1:dim
        proj[i] = dot(v, U[:, i])/dot(U[:,i], U[:,i])
    end
    return proj
end

function project_timeseries(data, U)
    flattened = reshape_data(data)
    timelen = size(data)[3]
    out = Array{Float64}(undef, (size(U)[2], timelen))
    for i in 1:timelen
        out[:,i] = projection(flattened[:,i], U) 
    end
    return out
end

# function back_to_data(v, U, S, V)
#     # v - vector of projection
#     # U, S, V - not necessarily full size, just as many columns as v is long 
#     v = normalize(v)
#     d = size(v)[1]
#     dims = (size(U)[1], size(V)[1], 1)
#     out = Array{Float64}(undef, dims)
#     # recover_mode = (U, S, V, i) ->  S[i] * reshape(U[:, i], (dims[1], 1)) * reshape(V[:, i], (1, dims[2])) 

#     for i in 1:d
#         # out = cat(out, recover_mode(U, S, V, i).*v[i], dims=3)
#         out = cat(out, U[:, i].*v[i], dims=3)
#     end
#     return abs.(sum(out, dims=3))
# end

function back_to_data(projts, basisU)
    return basisU*projts
end

function weighted_avg(ts::ncData)
    M, N, L = size(ts.data)
    Δϕ = reshape(2π / M * ones(M), (M, 1, 1))
    Δθ = reshape(π / N * ones(N) .* cos.(deg2rad.(ts.latvec)), (1, N, 1))
    # Δθ = reshape(π / N * ones(N) .* cos.(range(-π/2, π/2, N+2))[2:end-1], (1, N, 1))
    metric = (Δθ .* Δϕ) / (4π)
    return sum(metric .* ts.data, dims = (1,2))[:]
end
