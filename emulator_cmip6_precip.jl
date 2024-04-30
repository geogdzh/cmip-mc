# this is just an analysis script
using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list and the latvec to be used later on
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

##
scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]

# begin
#     monthtime = 1
#     fig = Figure(resolution=(2000,800))
#     ax = Axis(fig[1,1], title="true model variance")
#     heatmap!(ax, true_var[:,:,monthtime], colorrange=ext)
#     ax = Axis(fig[1,2], title="emulator variance")
#     heatmap!(ax, ens_vars_100[:,:,monthtime], colorrange=ext)
#     Colorbar(fig[1,3], colorrange=ext)
#     display(fig)
#     # save("figs/ens_var_compare_rcp85_month1.png", fig)
# end



variable = "temp"
numbers = [20, 100]
for scenario in scenarios[2:end]
    hfile = h5open("data/temp_precip/vars_$(variable)_$(scenario)_50ens.hdf5", "r")
    true_var = read(hfile, "true_var")
    true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
    close(hfile)

    wfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "cw") #let's make this also include means
    hfile = h5open("data/temp_precip/ens_vars_withpr_$(scenario).hdf5", "r") #this includes means
    for number in numbers
        ens_means = variable == "temp" ? read(hfile, "ens_means_tas_$(number)") : read(hfile, "ens_means_$(variable)_$(number)")
        ens_vars = variable == "temp" ? read(hfile, "ens_vars_tas_$(number)") : read(hfile, "ens_vars_$(variable)_$(number)")

        rmse = sqrt.(sum((true_var.-ens_vars).^2, dims=3)[:,:,1]./size(true_var)[3]) #time average rmse (shaped as a map)
        rmse_time = sqrt.(weighted_avg((true_var.-ens_vars).^2, latvec)) #spatial average rmse (shaped as a timeseries)
        rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
        rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries

        write(wfile, "rmse_vars_$(variable)_$(number)", rmse)
        write(wfile, "rmse_vars_time_$(variable)_$(number)", rmse_time)
        write(wfile, "rmse_means_$(variable)_$(number)", rmse_means)
        write(wfile, "rmse_means_time_$(variable)_$(number)", rmse_means_time)
    end
    close(hfile)
    close(wfile)
end


begin
    fig = Figure(resolution=(2000,1000))
    ax = Axis(fig[1,1], title="20 modes")
    # heatmap!(ax, pr_rmse_20, colorrange=pr_ext_var_rmse)
    heatmap!(ax, pr_rmse_means_20, colorrange=pr_ext_mean_rmse)
    ax = Axis(fig[1,2], title="100 modes")
    # heatmap!(ax, pr_rmse_100, colorrange=pr_ext_var_rmse)
    heatmap!(ax, pr_rmse_means_100, colorrange=pr_ext_mean_rmse)
    Colorbar(fig[1,3], colorrange=pr_ext_mean_rmse)
    display(fig)
    # save("figs/ens_mean_rmse_pr_ssp585.png", fig)
end

variable = "pr"
measure = "mean"
begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1], title="average RMSE of the mean", xlabel="years since 2014", ylabel="RMSE")
    linestyles = [ :dot, :solid]
    colors = [:red, :orange, :green, :blue, :purple]
    for (j, scenario) in enumerate(scenarios[2:end])
        if !(j in [2, 4])
            for (i, number) in enumerate(numbers)
                hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "r")
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_$(number)")
                lines!(ax, month_to_year_avg(ser), color=colors[j], alpha=0.6, linestyle=linestyles[i])
                close(hfile)
            end
        end
    end

    elems = [LineElement(color=:black, linestyle=:dot), LineElement(color=:black, linestyle=:solid)]
    elems_2 = [LineElement(color=:red), LineElement(color=:orange), LineElement(color=:green), LineElement(color=:blue), LineElement(color=:violet)]
    labels = ["20 modes",  "100 modes", "SSP119", "SSP126", "SSP245", "SSP370", "SSP585"]
    axislegend(ax, [elems..., elems_2...], labels, position=:lt)
    display(fig)
    # save("figs/ens_mean_rmse_time_pr_ssp_comparison_simplified.png", fig)
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
