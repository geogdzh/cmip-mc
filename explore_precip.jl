using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")

# goal: check the CMIP5 precip against ERA5 (but also just check ERA5 itself)

file_head = "/net/fs06/d3/ERA5/tprecip" #gonna use the monthly averaged big file
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/precip"

