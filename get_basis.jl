using LinearAlgebra, HDF5
include("utils.jl")

#### get basis - skip if loading it 

#first ensemble member of historical run
# file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_head = "data/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
ts = ncData(file_head*file_tail, "ts")
M, N, L1 = size(ts.data)
X = reshape_data(ts.data)

phfile = "/net/fs06/d3/CMIP5/MPI-GE/historical/precip/pr_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
pr = ncData(phfile, "pr")
Xp = reshape_data(pr.data)

#first ens member of 
# file_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/"
file_head = "data/RCP45/ts/"
file_tail = "ts_Amon_MPI-ESM_rcp45_r$(string(1, pad=3))i2005p3_200601-209912.nc"
ts2 = ncData(file_head*file_tail, "ts")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

phfile2 = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
pr2 = ncData(phfile2, "pr")
Xp2 = reshape_data(pr2.data)

fullX = hcat(X, X2)
U, S, V = svd(fullX)
d = 20
basis = U[:,1:d]

# hfile = h5open("data/temp_basis_20d.hdf5", "w")
# write(hfile, "basis", basis)
# close(hfile)

# fullX = hcat(vcat(X, Xp), vcat(X2, Xp2))
# U, S, V = svd(fullX)
# d = 100
# basis = U[:,1:d]
# hfile = h5open("data/temp-precip_basis.hdf5", "w")
# write(hfile, "basis", basis)
# close(hfile)
