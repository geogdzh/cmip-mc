using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
# include("simulator.jl")

# file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
# file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
# file = file_head*file_tail
# d = ncData(file, "ts")

# # ERA5 starts in 1959 ; historical ends in 2005.
# d_cut = d(1959,2005)
# reference_timevec = d_cut.timevec

##########
function get_normalization(file, var)
    if var == "t1000"
        var = "t"
    end
    scale_factor = ncgetatt(file, var, "scale_factor")
    add_offset = ncgetatt(file, var, "add_offset")
    return scale_factor, add_offset
end

function normalize_era5(x, scale_factor, add_offset)
    return (x * scale_factor) .+ add_offset
end

function get_month(date)
    return parse(Int64, Dates.format(date, "mm"))
end
function get_tvec(filename)
    return DateTime(1900,1,1)+Hour.(ncread(filename,"time"))
end
##########

erafile = "/net/fs06/d3/ERA5/pressure-levels/t1000-1959.nc"
sample_snap = ncread(erafile, "t", start=[1,1,1], count=[-1, -1, 1])
dims = size(dropdims(sample_snap, dims=3))
goal_data_array = Array{Float64}(undef, dims)
#cycle through each year and average over a month:
begin
    current_date = Dates.DateTime(1959,1,1)
    working_avg = Array{Float64}(undef, dims)
    for yr in 1959:1960#2005
        println("starting year "*string(yr))
        file = "/net/fs06/d3/ERA5/pressure-levels/t1000-$yr.nc" 
        tvec = get_tvec(file)
        scale_factor, add_offset = get_normalization(file, "t")
        for ind in eachindex(tvec)
            date = tvec[ind]
            month = get_month(date)
            snap = dropdims(ncread(file, "t", start=[1,1,ind], count=[-1, -1, 1]), dims=3)
            snap = normalize_era5(snap, scale_factor, add_offset)
            if month == get_month(current_date)
                working_avg = cat(working_avg, snap, dims=3)
            else
                println("starting month "*string(month))
                # run the averaging
                avg = mean(working_avg, dims=3)
                # add it to the right place in ref vec
                goal_data_array = cat(goal_data_array, avg, dims=3)
                # reset working average
                working_avg = snap 
            end
            current_date = date
        end

    end
end

using JLD

JLD.save("data/test_era5_2.jld", "goal_data_array", goal_data_array)