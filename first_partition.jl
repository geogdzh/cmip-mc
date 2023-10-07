using NetCDF
using Dates
using Statistics
using CairoMakie

include("utils.jl")

# select reference values

file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"

refs = Dict(i => [] for i in 1:12)
for i in [3, 47, 86]
    file_tail = "ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc"
    file = file_head*file_tail
    d = ncData(file, "ts")#(-20, 15, 5, 33) #CHANGE TO CHANGE SLICE
    for month in 1:12  
        push!(refs[month], d.data[:,:,12])    #worth noting this makes it from just the first year
    end
end

heatmap(refs[1][1])
heatmap(refs[1][2])
heatmap(refs[1][3])

# now we have a monthly dict of three reference images
# so, next, want to "walk through each simulation and get a markov chain (one for each month)


include("simulator.jl")

mc_dict = Dict(i => [] for i in 1:12)

for i in 1:3#100
    file_tail = "ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc"
    file = file_head*file_tail
    d = ncData(file, "ts")
    mcs = step_through(d, refs)
    for j in 1:12
        push!(mc_dict[j], mcs[j])
    end
end

# we're left with, for each month, 100 mc's from each ensemble member.
# want to get an operator: generator or PF ? (using alll of the data to also get an uncertainty estimate)
# let's start with generator because we can

