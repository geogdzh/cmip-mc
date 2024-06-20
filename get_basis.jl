using LinearAlgebra, HDF5
include("utils.jl")
include("eof_util.jl")

#### get basis - based ONLY on one ens memeber (should this be substantiated?)
non_dim = false
use_metrics = true

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
# file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/" ### CHANGE HERE FOR DIFF LOCATION
parent_folder = non_dim ? "nondim" : "temp_precip"
parent_folder = use_metrics ? "metrics" : parent_folder

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

phfile = file_head*"historical/pr/r10i1p1f1_historical_pr.nc"
pr = ncData(phfile, "pr")
Xp = reshape_data(use_metrics ? sqrt.(metric) .* pr.data : pr.data)


#first ens member of the model run
file_tail = "ssp585/tas/r10i1p1f1_ssp585_tas.nc"
ts2 = ncData(file_head*file_tail, "tas")
X2 = reshape_data(ts2.data)
M, N, L2 = size(use_metrics ? sqrt.(metric) .* ts2.data : ts2.data)

phfile2 = file_head*"ssp585/pr/r10i1p1f1_ssp585_pr.nc"
pr2 = ncData(phfile2, "pr")
Xp2 = reshape_data(use_metrics ? sqrt.(metric) .* pr2.data : pr2.data)

if non_dim
    temp_factor = maximum(X)
    X = X ./ temp_factor
    X2 = X2 ./ temp_factor
    pr_factor = maximum(Xp)
    Xp = Xp ./ pr_factor
    Xp2 = Xp2 ./ pr_factor
end

fullX = hcat(X, X2)
fullXp = hcat(Xp, Xp2)
full = vcat(fullX, fullXp)
U, S, V = svd(full)
d = 1000
basis = U[:,1:d]

hfile = h5open("data/$(parent_folder)/temp_precip_basis_$(d)d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
if non_dim
    write(hfile, "temp_factor", temp_factor)
    write(hfile, "pr_factor", pr_factor)
end
if use_metrics
    write(hfile, "metric", metric)
end
close(hfile)