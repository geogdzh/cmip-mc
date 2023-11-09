
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

# function get_monthly_eofs(d::ncData)
#     refs = Dict(i => [] for i in 1:12)
#     for month in 1:12
#         X = reshape_data(d.data[:,:,month:12:end])
#         U, S, V = svd(X) #columns of U span the column space of X
#         refs[month] = U
#     end
#     return refs
# end

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

function back_to_data(v, U, S, V)
    # v - vector of projection
    # U, S, V - not necessarily full size, just as many columns as v is long 
    v = normalize(v)
    d = size(v)[1]
    dims = (size(U)[1], size(V)[1], 1)
    out = Array{Float64}(undef, dims)
    # recover_mode = (U, S, V, i) ->  S[i] * reshape(U[:, i], (dims[1], 1)) * reshape(V[:, i], (1, dims[2])) 

    for i in 1:d
        # out = cat(out, recover_mode(U, S, V, i).*v[i], dims=3)
        out = cat(out, U[:, i].*v[i], dims=3)
    end
    return abs.(sum(out, dims=3))
end
