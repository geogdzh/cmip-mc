using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
include("simulator.jl")

all_eofs_cat = JLD.load("data/all_eofs_historical.jld", "all_eofs_cat")
# all_eofs = [all_eofs_cat[:,:,i] for i in 1:100]

M = 192
N = 96

U1 = all_eofs_cat[:,:,1] #U for first ensemble memeber

eof1 = shape_data(U1[:,1], M, N)

heatmap(eof1)
extrema(eof1)

basis = U1[:,1:5]

#sample snapshot:
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail
ts = ncData(file, "ts")
v = ts.data[:,:, 1]

# projection(v, basis)