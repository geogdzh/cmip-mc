
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
    dim = size(U,2)
    proj = zeros(dim)
    for i in 1:dim
        proj[i] = dot(v, U[:, i])/dot(U[:,i], U[:,i])
    end
    return proj
end

function project_timeseries()
    #pass
end