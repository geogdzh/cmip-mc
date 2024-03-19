using CairoMakie, ProgressBars, HDF5, NetCDF, NCDatasets
include("utils.jl")
include("eof_util.jl")

# goal: check the CMIP5 precip against ERA5 (but also just check ERA5 itself)

era5_file_head = "/net/fs06/d3/ERA5/" #gonna use the monthly averaged big file
cmip5_file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/precip/"

era5 = h5open(era5_file_head*"monthly_averaged_tprecip/era5_monthly_avg_tprecip_1959_2002.h5", "r")
data_avg = read(era5, "monthly_averaged_precip")
timevec = read(era5, "time")
close(era5)

# era5 = h5open(era5_file_head*"monthly_averaged_t2m/era5_monthly_avg_1960_2005.h5", "r")
# data = read(era5, "monthly_averaged_temperature")
# timevec = read(era5, "time")
# close(era5)

# for i in eachindex(timevec)
#     if mod(i, 12) in [1,3,5,7,8,10,12]
#         data[:,:,i] .*= 24*31
#     elseif mod(i, 12) in [4,6,9,11]
#         data[:,:,i] .*= 24*30
#     else
#         data[:,:,i] .*= 24*28
#     end
# end

era5 = h5open("era5_monthly_total_tprecip_1959_2002.h5", "r")
data_total = read(era5, "monthly_total_precip")
timevec = read(era5, "time")
close(era5)

sample_file = era5_file_head*"tprecip/tprecip-1959.nc"
lonvec = ncread(sample_file, "longitude")
latvec = ncread(sample_file, "latitude")
avg_era5 = ncData(data_avg, lonvec, latvec, timevec)
total_era5 = ncData(data_total, lonvec, latvec, timevec)

avg = weighted_avg(total_era5) #assuming this returns a timeseries

begin
    fig = Figure(resolution=(1200, 800))
    ax = Axis(fig[1, 1], xlabel="Year", ylabel="Precipitation", title="ERA5 Precipitation")
    lines!(ax, avg)
    display(fig)
end

begin
    fig = Figure(resolution=(1200, 800))
    ax = Axis(fig[1,1])
    heatmap!(ax, monthly_era5.data[:,:,1])
    display(fig)
end

# ok let's start over with just the original ERA5 data
file = era5_file_head*"tprecip/tprecip-1959.nc"
time = ncread(file, "time")

ds = Dataset(file)
var = "tp"
# maximum(ds[var][:, :, 1].*scale_factor .+ add_offset)
tmp = zeros(1440, 721, 24)
for i in 1:24
    tmp[:, :, 1] += ds[var][:, :, i]
end
heatmap(tmp[:,:,1])

scale_factor = ncgetatt(file, var, "scale_factor")
add_offset = ncgetatt(file, var, "add_offset")

###


