# actual tool used to generate era5 monthly averages

using NCDatasets, ProgressBars
# data_directory = "/net/fs06/d3/ERA5/t2m"
data_directory = "/net/fs06/d3/ERA5/tprecip"

years = [x for x in 1959:2002]
monthly_averaged_var = zeros(1440, 721, 12*length(years));
for j in eachindex(years)
    println("starting year "*string(years[j]))
    # file = "t2m-" * string(years[j]) * ".nc"
    file = "tprecip-" * string(years[j]) * ".nc"
    filepath = joinpath(data_directory, file)
    ds = Dataset(filepath)
    # var = "t2m"
    var = "tp"

    running_index = zeros(12); #temp
    for i in ProgressBar(eachindex(ds["time"]))
        tval = ds[var][:, :, i]
        month_index = NCDatasets.Dates.month(ds["time"][i])
        running_index[month_index] += 1
        monthly_averaged_var[:, :, 12*(j-1)+month_index] .+= tval
    end
    for i in 1:12
        monthly_averaged_var[:, :, 12*(j-1)+i] .*= 1/running_index[i]
    end
end

timevec = Vector{String}()
for year in years
    push!(timevec, [string(year)*"-"*(string(i, pad=2)) for i in 1:12]...)
end

using HDF5
h5file = h5open("era5_monthly_avg_tprecip_1960_2005.h5", "w") 
write(h5file, "monthly_averaged_precip", monthly_averaged_var)
write(h5file, "time", timevec)
close(h5file)