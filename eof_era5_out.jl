using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, HDF5
include("utils.jl")
include("eof_util.jl")


era5 = h5open("/net/fs06/d3/ERA5/monthly_averaged_t2m/era5_monthly_avg_1960_2005.h5", "r")
era5_temp = read(era5, "monthly_averaged_temperature")
time = read(era5, "time")
close(era5)
Me, Ne, Le = size(era5_avg_temp)

function average_data(temperature)
    M, N, L = size(temperature)
    pass1 = zeros(180, N, L)
    for i in 1:180
        pass1[i, :, :] .= mean(temperature[1 + (i-1)*8:i * 8, :, :], dims = 1)[1, :, :]
    end
    pass2 = zeros(180, 90, L)
    for i in 1:90
        pass2[:, i, :] .= mean(pass1[:, 1 + (i-1)*8:i * 8, :], dims = 2)[:,1, :]
    end
    return pass2
end

era5_avg_temp = average_data(era5_temp)


X = reshape_data(era5_avg_temp)[:,1+30*12:end]
U, S, V = svd(X)
# heatmap(shape_data(U[:,1], Me, Ne)[:,end:-1:1], reverse=true)

file = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
cmip = ncData(file, "ts")(1960,2005)
cmip_temp = cmip.data
Mc, Nc, Lc = size(cmip_temp)
Xc = reshape_data(cmip_temp)[:,1+30*12:end]
Uc, Sc, Vc = svd(Xc)

# all_ext = []

begin
    fig = Figure(resolution=(1200, 2000))
    for i in 1:5
        mode = shape_data(U[:,i], Me, Ne)[:,end:-1:1]
        # ext = extrema(mode)
        # push!(all_ext, ext)
        ax1 = Axis(fig[i, 1], title= i==1 ? "ERA5" : "")
        heatmap!(ax1, mode, colorrange=all_ext[i])
        mode2 = shape_data(Uc[:,i],Mc, Nc)
        ax2 = Axis(fig[i, 2], title= i==1 ? "MPI-GE (1 ensemble member)" : "")
        heatmap!(ax2, mode2, colorrange=all_ext[i])
    end
    # title(fig, "ERA5 vs. MPI-GE modes")
    display(fig)
    save("figs/mode_comparison_1990-2005.png", fig)
end