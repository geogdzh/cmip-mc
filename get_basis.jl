using LinearAlgebra, HDF5
include("utils.jl")
include("eof_util.jl")

#### get basis - based ONLY on one ens memeber (should this be substantiated?)

#first ensemble member of historical run

# file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/" ### CHANGE HERE FOR DIFF LOCATION

file_tail = "historical/tas/r10i1p1f1_historical_tas.nc"
ts = ncData(file_head*file_tail, "tas")
M, N, L1 = size(ts.data)
X = reshape_data(ts.data)

phfile = file_head*"historical/pr/r10i1p1f1_historical_pr.nc"
pr = ncData(phfile, "pr")
Xp = reshape_data(pr.data)

#first ens member of the model run
file_tail = "ssp585/tas/r10i1p1f1_ssp585_tas.nc"
ts2 = ncData(file_head*file_tail, "tas")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

phfile2 = file_head*"ssp585/pr/r10i1p1f1_ssp585_pr.nc"
pr2 = ncData(phfile2, "pr")
Xp2 = reshape_data(pr2.data)

fullX = hcat(X, X2)
fullXp = hcat(Xp, Xp2)
full = vcat(fullX, fullXp)
U, S, V = svd(full)
d = 1000
basis = U[:,1:d]

hfile = h5open("data/temp_precip/temp_precip_basis_1000d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
close(hfile)