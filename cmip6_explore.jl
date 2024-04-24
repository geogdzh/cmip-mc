using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 


scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]
scenario = scenarios[1]
ens_ind = 8
for ens_ind in 1:50
    println(ens_ind)
    filehead = "/net/fs06/d3/mgeo/CMIP6/interim/"
    filetail = "$(scenario)/tas/r$(ens_ind)i1p1f1_$(scenario)_tas.nc"
    file = filehead * filetail

    tas = ncData(file, "tas")
    # tas.timevec
    # timeseries = month_to_year_avg(weighted_avg(tas))

    println(count_nan(tas.data)  )
end

lines(timeseries)
heatmap(tas.data[:,:,2])

function count_nan(arr::AbstractArray)
    count = 0
    locs = []
    for i in 1:M
        for j in 1:N
            for k in 1:L1
                value = arr[i,j,k]
                if ismissing(value)
                    count += 1
                    push!(locs, (i,j,k))
                end
            end
        end
    end
    return count, locs
end

function count_nan(arr::AbstractArray)
    count = 0
    for value in arr
        if ismissing(value)
            count += 1
        end
    end
    return count
end

cnt, locs = count_nan(tas.data) #221184


hfile = h5open("data/CMIP5/training_data_rcp85_20d_90ens.hdf5", "r")
projts_20 = read(hfile, "projts")
ens_gmt_20 = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/CMIP5/training_data_rcp85_100d_90ens.hdf5", "r")
projts_100 = read(hfile, "projts")
ens_gmt_100 = read(hfile, "ens_gmt")
close(hfile)

ens_gmt_20 == ens_gmt_100
projts_100[1:20,:,:] == projts_20

scenario = "ssp585"
ds = Dataset("/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/tas/r$(1)i1p1f1_$(scenario)_tas.nc")
tas = ds["tas"][:]
close(ds)

tas = ncData("/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/tas/r$(10)i1p1f1_$(scenario)_tas.nc", "tas")
gmt_119 = get_gmt_list(tas)
gmt_585 = get_gmt_list(tas)
lines(gmt_119)
lines(gmt_585)

#

file = "/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/r1i1p1f1/historical/pr/250_km/mon/1850/CMIP6_MPI-ESM1-2-LR_r1i1p1f1_historical_pr_250_km_mon_gn_1850.nc"
pr = ncData(file, "pr")

using NCDatasets
ds = Dataset(file)

ncgen(file, "samplepr.jl")