using LinearAlgebra, HDF5
include("utils.jl")

#### get basis - based ONLY on one ens memeber (should this be substantiated?)

#first ensemble member of historical run
# file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
# file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/historical/tas/"
file_tail = "r10i1p1f1_historical_tas.nc"
ts = ncData(file_head*file_tail, "tas")
M, N, L1 = size(ts.data)
X = reshape_data(ts.data)

# phfile = "/net/fs06/d3/CMIP5/MPI-GE/historical/precip/pr_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
phfile = "/net/fs06/d3/mgeo/CMIP6/interim/historical/pr/r10i1p1f1_historical_pr.nc"
pr = ncData(phfile, "pr")
Xp = reshape_data(pr.data)

#first ens member of the model run
# file_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP85/ts/"
# file_tail = "ts_Amon_MPI-ESM_rcp85_r$(string(1, pad=3))i2005p3_200601-209912.nc"
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/ssp585/tas/"
file_tail = "r10i1p1f1_ssp585_tas.nc"
ts2 = ncData(file_head*file_tail, "tas")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

# phfile2 = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
phfile2 = "/net/fs06/d3/mgeo/CMIP6/interim/ssp585/pr/r10i1p1f1_ssp585_pr.nc"
pr2 = ncData(phfile2, "pr")
Xp2 = reshape_data(pr2.data)

fullX = hcat(X, X2)
fullXp = hcat(Xp, Xp2)
full = vcat(fullX, fullXp)
U, S, V = svd(full)
d = 1000
basis = U[:,1:d]

hfile = h5open("data/temp_precip_basis_1000d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
close(hfile)

# fullX = hcat(vcat(X, Xp), vcat(X2, Xp2))
# U, S, V = svd(fullX)
# d = 100
# basis = U[:,1:d]
# hfile = h5open("data/temp-precip_basis.hdf5", "w")
# write(hfile, "basis", basis)
# close(hfile)
