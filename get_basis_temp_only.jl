using LinearAlgebra, HDF5
include("utils.jl")
include("eof_util.jl")

#### get basis - based ONLY on one ens memeber (should this be substantiated?)
non_dim = false
use_metrics = false

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
# file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/" ### CHANGE HERE FOR DIFF LOCATION
parent_folder = "temp_only"

#first ensemble member of historical run
file_tail = "historical/tas/r10i1p1f1_historical_tas.nc"
ts = ncData(file_head*file_tail, "tas")
M, N, L1 = size(ts.data)

if use_metrics
    latvec = ts.latvec
    Δϕ = reshape(2π / M * ones(M), (M, 1, 1))
    Δθ = reshape(π / N * ones(N) .* cos.(deg2rad.(latvec)), (1, N, 1))
    metric = (Δθ .* Δϕ) / (4π)
end

X = reshape_data(use_metrics ? sqrt.(metric) .* ts.data : ts.data)

#first ens member of the model run
file_tail = "ssp585/tas/r10i1p1f1_ssp585_tas.nc"
ts2 = ncData(file_head*file_tail, "tas")
X2 = reshape_data(ts2.data)
M, N, L2 = size(use_metrics ? sqrt.(metric) .* ts2.data : ts2.data)

if non_dim
    temp_factor = maximum(X)
    X = X ./ temp_factor
    X2 = X2 ./ temp_factor
end

full = hcat(X, X2)
U, S, V = svd(full)
d = 2000
basis = U[:,1:d]

hfile = h5open("data/$(parent_folder)/temp_basis_$(d)d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
write(hfile, "S", S)
if non_dim
    write(hfile, "temp_factor", temp_factor)
end
if use_metrics
    write(hfile, "metric", metric)
end
close(hfile)