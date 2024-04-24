# this is just an analysis script
using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list
scenario = "ssp585"
file3 = "/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/tas/r$(1)i1p1f1_$(scenario)_tas.nc"
ts3 = ncData(file3, "tas") 
gmt_list = get_gmt_list(ts3)
M, N, L = size(ts3.data)
latvec = ts3.latvec


# sim = emulate(gmt_list, mean_coefs, chol_coefs)
# sim = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
# newdata = back_to_data(sim, basis)
# simts = ncData(shape_data(newdata, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
# sim_gmt = get_gmt_list(simts)

#split the data
# newts = newdata[1:M*N,:]
# newtp = newdata[M*N+1:end, :]

#### correct way of comparing rmse -- only needs chol_coefs as the 'emulator' component! 

hfile = h5open("data/vars_ssp585_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
true_var = read(hfile, "true_var")
num_ens_members = read(hfile, "num_ens_members")
true_ens_gmt = read(hfile, "true_ens_gmt")
true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
close(hfile)
# ext = (0., extrema(true_var)[2])


### emulator vars adn means

function get_ens_vars(d, true_ens_gmt; get_means=false) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = h5open("data/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    close(hfile)
    ens_vars = zeros(M, N, L)
    for m in ProgressBar(1:Int(L/12))
        for n in 1:12
            if get_means
                ens_vars[:,:,(m-1)*12+n] = shape_data(back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1], basis), M, N)
            else
                co = get_cov(true_ens_gmt[m], chol_coefs)[:,:,n] #there is maybe a more efficient way to do this?
                ens_vars[:,:,(m-1)*12+n] = shape_data(sum([co[i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]), M, N)
            end
        end
    end
    return ens_vars
end

######################################################################################################################
# ens_means_20 = get_ens_vars(20, true_ens_gmt; get_means=true)
# ens_means_60 = get_ens_vars(60, true_ens_gmt; get_means=true)
# ens_means_100 = get_ens_vars(100, true_ens_gmt; get_means=true)
# ens_means_200 = get_ens_vars(200, true_ens_gmt; get_means=true)

# rmse_means_20 = sqrt.(sum((true_ens_mean.-ens_means_20).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map) #this is correct
# # ext_means = (0., maximum(rmse_means_20)) #this is FIXED based on teh ssp585 compariosn
# rmse_means_60 = sqrt.(sum((true_ens_mean.-ens_means_60).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])
# rmse_means_100 = sqrt.(sum((true_ens_mean.-ens_means_100).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])
# rmse_means_200 = sqrt.(sum((true_ens_mean.-ens_means_200).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])

# rmse_means_time_20 = sqrt.(weighted_avg((true_ens_mean.-ens_means_20).^2, latvec)) #spatial average rmse (shaped as a timeseries)
# rmse_means_time_60 = sqrt.(weighted_avg((true_ens_mean.-ens_means_60).^2, latvec))
# rmse_means_time_100 = sqrt.(weighted_avg((true_ens_mean.-ens_means_100).^2, latvec))
# rmse_means_time_200 = sqrt.(weighted_avg((true_ens_mean.-ens_means_200).^2, latvec))

# rmse_mean_numbers = zeros((4, 5))
# all_rmse_means_time = zeros((4, 5, L))
# hfile = h5open("data/ens_means_rmse_comparison.hdf5", "w")
# write(hfile, "rmse_mean_numbers", rmse_mean_numbers)
# write(hfile, "all_rmse_means_time", all_rmse_means_time)
# close(hfile)

# ssp119_time_avgs = (rmse_means_time_20, rmse_means_time_60, rmse_means_time_100, rmse_means_time_200)
# ssp126_time_avgs = (rmse_means_time_20, rmse_means_time_60, rmse_means_time_100, rmse_means_time_200)
# ssp245_time_avgs = (rmse_means_time_20, rmse_means_time_60, rmse_means_time_100, rmse_means_time_200)
# ssp370_time_avgs = (rmse_means_time_20, rmse_means_time_60, rmse_means_time_100, rmse_means_time_200)
# ssp585_time_avgs = (rmse_means_time_20, rmse_means_time_60, rmse_means_time_100, rmse_means_time_200)
# all_time_avgs = (ssp119_time_avgs, ssp126_time_avgs, ssp245_time_avgs, ssp370_time_avgs, ssp585_time_avgs)
######################################################################################################################

# ens_vars_20 = get_ens_vars(20, true_ens_gmt)
# ens_vars_60 = get_ens_vars(60, true_ens_gmt)
# ens_vars_100 = get_ens_vars(100, true_ens_gmt)
# ens_vars_200 = get_ens_vars(200, true_ens_gmt)

# hfile = h5open("data/ens_vars_ssp585.hdf5", "w")
# write(hfile, "ens_vars_20", ens_vars_20)
# write(hfile, "ens_vars_60", ens_vars_60)
# write(hfile, "ens_vars_100", ens_vars_100)
# write(hfile, "ens_vars_200", ens_vars_200)
# close(hfile)

##
scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]
scenario = scenarios[2]

hfile = h5open("data/vars_$(scenario)_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
true_var = read(hfile, "true_var")
num_ens_members = read(hfile, "num_ens_members")
true_ens_gmt = read(hfile, "true_ens_gmt")
true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
close(hfile)

hfile = h5open("data/ens_vars_$(scenario).hdf5", "r")
ens_vars_20 = read(hfile, "ens_vars_20")
ens_vars_60 = read(hfile, "ens_vars_60d")
ens_vars_100 = read(hfile, "ens_vars_100d")
ens_vars_200 = read(hfile, "ens_vars_200d")
close(hfile)

begin
    monthtime = 1
    fig = Figure(resolution=(2000,800))
    ax = Axis(fig[1,1], title="true model variance")
    heatmap!(ax, true_var[:,:,monthtime], colorrange=ext)
    ax = Axis(fig[1,2], title="emulator variance")
    heatmap!(ax, ens_vars_100[:,:,monthtime], colorrange=ext)
    Colorbar(fig[1,3], colorrange=ext)
    display(fig)
    # save("figs/ens_var_compare_rcp85_month1.png", fig)
end

#to be clear, right now we're looking at the RMSE of the variance! i.e. how well does the emulator capture model variance

rmse_20 = sqrt.(sum((true_var.-ens_vars_20).^2, dims=3)[:,:,1]./size(true_var)[3]) #time average rmse (shaped as a map)
# ext_var_rmse = (0., maximum(rmse_20)) # - fix this based on the ssp585 comparison
rmse_60 = sqrt.(sum((true_var.-ens_vars_60).^2, dims=3)[:,:,1]./size(true_var)[3])
rmse_100 = sqrt.(sum((true_var.-ens_vars_100).^2, dims=3)[:,:,1]./size(true_var)[3])
rmse_200 = sqrt.(sum((true_var.-ens_vars_200).^2, dims=3)[:,:,1]./size(true_var)[3])

rmse_time_20 = sqrt.(weighted_avg((true_var.-ens_vars_20).^2, latvec)) #spatial average rmse (shaped as a timeseries
rmse_time_60 = sqrt.(weighted_avg((true_var.-ens_vars_60).^2, latvec))
rmse_time_100 = sqrt.(weighted_avg((true_var.-ens_vars_100).^2, latvec))
rmse_time_200 = sqrt.(weighted_avg((true_var.-ens_vars_200).^2, latvec))

# ssp119_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# ssp126_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# ssp245_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# ssp370_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# ssp585_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# all_time_avgs = (ssp119_time_avgs, ssp126_time_avgs, ssp245_time_avgs, ssp370_time_avgs, ssp585_time_avgs)

# rmse_var_numbers = zeros((4, 5))
# all_rmse_vars_time = zeros((4, 5, L))
# hfile = h5open("data/ens_vars_rmse_comparison.hdf5", "w")
# write(hfile, "rmse_var_numbers", rmse_var_numbers)
# write(hfile, "all_rmse_vars_time", all_rmse_vars_time)
# close(hfile)


begin
    fig = Figure(resolution=(2000,1500))
    ax = Axis(fig[1,1], title="20 modes")
    heatmap!(ax, rmse_20, colorrange=ext_var_rmse)
    # heatmap!(ax, rmse_means_20, colorrange=ext_means)
    ax = Axis(fig[1,2], title="60 modes")
    heatmap!(ax, rmse_60,  colorrange=ext_var_rmse)
    # heatmap!(ax, rmse_means_60,  colorrange=ext_means)
    ax = Axis(fig[2,1], title="100 modes")
    heatmap!(ax, rmse_100, colorrange=ext_var_rmse)
    # heatmap!(ax, rmse_means_100, colorrange=ext_means)
    ax = Axis(fig[2,2], title="200 modes")
    heatmap!(ax, rmse_200, colorrange=ext_var_rmse)
    # heatmap!(ax, rmse_means_200, colorrange=ext_means)
    Colorbar(fig[1,3], colorrange=ext_var_rmse)
    display(fig)
    # save("figs/ens_var_rmse_ssp585.png", fig)
end


begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1], title="average RMSE of the variance", xlabel="years since 2014", ylabel="RMSE")
    linestyles = [ :dot, :dashdot, :dash, :solid]
    colors = [:red, :orange, :green, :blue, :purple]
    for (j, scenario) in enumerate(all_time_avgs)
        if !(j in [2, 4])
            for (i, ser) in enumerate(scenario)
                # rmse_mean_numbers[i, j] = mean(ser)
                # all_rmse_means_time[i, j, :] = ser
                # rmse_var_numbers[i, j] = mean(ser)
                # all_rmse_vars_time[i, j, :] = ser
                lines!(ax, month_to_year_avg(ser), color=colors[j], alpha=0.6, linestyle=linestyles[i])
            end
        end
    end
    # for (i, ser) in enumerate(rcp26_time_avgs)
    #     rmse_var_numbers[i, 1] = mean(ser)
    #     all_rmse_vars_time[i, 1, :] = ser
    #     lines!(ax, month_to_year_avg(ser), color=:orange, alpha=0.6, linestyle=linestyles[i])
    # end
    # for (i, ser) in enumerate(rcp45_time_avgs)
    #     rmse_var_numbers[i, 2] = mean(ser)
    #     all_rmse_vars_time[i, 2, :] = ser
    #     lines!(ax, month_to_year_avg(ser), color=:green, alpha=0.6, linestyle=linestyles[i])
    # end
    # for (i, ser) in enumerate(rcp85_time_avgs)
    #     rmse_var_numbers[i, 3] = mean(ser)
    #     all_rmse_vars_time[i, 3, :] = ser
    #     lines!(ax, month_to_year_avg(ser), color=:blue, alpha=0.6, linestyle=linestyles[i])
    # end
    # lines!(ax, month_to_year_avg(rmse_means_time_20), label="20 modes")
    # lines!(ax, month_to_year_avg(rmse_means_time_60), label="60 modes")
    # lines!(ax, month_to_year_avg(rmse_means_time_100), label="100 modes")
    # lines!(ax, month_to_year_avg(rmse_means_time_200), label="200 modes")

    elems = [LineElement(color=:black, linestyle=:dot), LineElement(color=:black, linestyle=:dashdot), LineElement(color=:black, linestyle=:dash), LineElement(color=:black, linestyle=:solid)]
    elems_2 = [LineElement(color=:red), LineElement(color=:orange), LineElement(color=:green), LineElement(color=:blue), LineElement(color=:violet)]
    labels = ["20 modes", "60 modes", "100 modes", "200 modes", "SSP119", "SSP126", "SSP245", "SSP370", "SSP585"]
    axislegend(ax, [elems..., elems_2...], labels, position=:lb)
    display(fig)
    # save("figs/ens_var_rmse_time_ssp_comparison_simplified.png", fig)
end






####### compare gmt representation w some ensemble members

hfile = h5open("data/gaussian_emulator_ssp585_$(20)d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
close(hfile)

ens_mem = 10
gmts = zeros((ens_mem, Int(L/12)))
gmts_cov = zeros((ens_mem, Int(L/12)))
for i in ProgressBar(1:ens_mem)
    s = emulate(gmt_list, mean_coefs, chol_coefs)
    # s = emulate(gmt_list, mean_coefs, chol_coefs; no_cov=true)
    data = back_to_data(s, basis)[1:M*N, :]
    sts = ncData(shape_data(data, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
    gmts_cov[i,:] = get_gmt_list(sts)
    # gmts[i,:] = get_gmt_list(sts)
end
avg = dropdims(mean(gmts, dims=1),dims=1)
avg_cov = dropdims(mean(gmts_cov, dims=1),dims=1)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    # lines!(ax,  sim_gmt, label="emulator")
    for i in 1:ens_mem
        lines!(ax, gmts_cov[i, :], color=:blue, alpha=0.2) #blue = with covariance accounted for ?
        # lines!(ax, gmts[i, :], color=:orange, alpha=0.3) #orange = without covariance accounted for
    end
    # lines!(ax, avg, label="emulator ensemble avg", color=:red)
    lines!(ax, avg_cov, label="emulator ensemble avg with cov", color=:blue)
    lines!(ax,  gmt_list, label="model", color=:black)
    axislegend(ax, position=:lt)
    display(fig)
    # save("figs/sim_gmt_compare_cov_rcp85.png", fig)
end
