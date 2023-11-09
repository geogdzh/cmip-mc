using NCDatasets, ProgressBars
data_directory = "/net/fs06/d3/ERA5/t2m"

years = [x for x in 1960:2005]
monthly_averaged_temperature = zeros(1440, 721, 12*length(years));
for j in eachindex(years)
    println("starting year "*string(years[j]))
    file = "t2m-" * string(years[j]) * ".nc"
    filepath = joinpath(data_directory, file)
    ds = Dataset(filepath)
    var = "t2m"

    running_index = zeros(12); #temp
    for i in ProgressBar(eachindex(ds["time"]))
        tval = ds[var][:, :, i]
        month_index = NCDatasets.Dates.month(ds["time"][i])
        running_index[month_index] += 1
        monthly_averaged_temperature[:, :, 12*(j-1)+month_index] .+= tval
    end
    for i in 1:12
        monthly_averaged_temperature[:, :, 12*(j-1)+i] .*= 1/running_index[i]
    end
end

timevec = Vector{String}()
for year in years
    push!(timevec, [string(year)*"-"*(string(i, pad=2)) for i in 1:12]...)
end


# using JLD
# JLD.save("data/era5_monthly_avg_historical.jld", "monthly_averaged_temperature", monthly_averaged_temperature)

using HDF5
h5file = h5open("era5_monthly_avg_1960_2005.h5", "w") 
write(h5file, "monthly_averaged_temperature", monthly_averaged_temperature)
write(h5file, "time", timevec)
close(h5file)