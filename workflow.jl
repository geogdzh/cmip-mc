using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
include("simulator.jl")

#### Step 1: define reference

#first ensemble member of historical run
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail
ts = ncData(file, "ts")

#get EOF basis
X = reshape_data(ts.data)
U, S, V = svd(X)
basis = U[:,1:5]

#partition state space
