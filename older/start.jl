using NetCDF
using Dates
using Statistics
using CairoMakie

include("utils.jl")

file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r001i1850p3_185001-200512.nc" #number after r = ensemble memeber out of 100
file = file_head*file_tail

ncinfo(file)

ts = ncread(file, "ts")
slice = ts[:,:,1]

heatmap(slice)

timevec = get_tvec(file)

latvec = ncread(file, "lat")
lonvec = ncread(file, "lon")

nigeria = slice[find_closest(lonvec, 2):find_closest(lonvec, 15), find_closest(latvec, 4):find_closest(latvec, 14)]
tester = slice[find_closest(lonvec, 125):find_closest(lonvec, 151), find_closest(latvec, 30):find_closest(latvec, 46)]

heatmap(tester)


west_africa = slice_map(ts, lonvec, latvec, -20, 15, 0, 30)
heatmap(west_africa[:,:,1])

ts = ncData(file, "ts")

ts.latvec

# ts_slice = slice_map(ts, -20, 15, 0, 30)
ts_slice = ts(-85, -65, 35, 48) #Eastern US coast sorta
heatmap(ts_slice.data[:,:,1])

heatmap(ts.data[:,:,1])

month(ts.timevec[1])

