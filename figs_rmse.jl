using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list and the latvec to be used later on
file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

##
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 #for CMIP6
l1, l2 = 165, 86

variable = "temp"
numbers = [20, 100]
for scenario in scenarios[2:end]
    # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

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

numbers = [20, 100]
ks = [x for x in 1:3]
variable = "pr"
measure = "mean"
testing_k = false #iMPLEMENT THIS HERE
begin
    fig = Figure(resolution=(800,600))
    ax = Axis(fig[1,1], title="Average RMSE of the ensemble $(measure) for $(variable == "temp" ? "temperature" : "precipitation")", xlabel="Year", ylabel="RMSE")
    linestyles = [ :dot, :solid]
    for (j, scenario) in enumerate(scenarios[2:end])
        for (i, number) in enumerate(numbers)
            hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "r")
            ser = read(hfile, "rmse_$(measure)s_time_$(variable)_$(number)")
            lines!(ax,time_future, month_to_year_avg(ser), color=scenario_colors[scenario], alpha=0.6, linestyle=linestyles[i])
            close(hfile)
        end
    end

    elems = [LineElement(color=:black, linestyle=:dot), LineElement(color=:black, linestyle=:solid)]
    elems_2 = [ LineElement(color=scenario_colors["ssp119"]), LineElement(color=scenario_colors["ssp245"]), LineElement(color=scenario_colors["ssp585"])]
    labels = ["20 modes",  "100 modes", "SSP119", "SSP245", "SSP585"]
    axislegend(ax, [elems..., elems_2...], labels, position=:lt)
    display(fig)
    # save("figs/ens_$(measure)_rmse_time_$(variable)_ssp_comparison.png", fig)
end

variable = "pr"
measure = "var"
begin
    fig = Figure(resolution=(800,600))
    ax = GeoAxis(fig[1,1], title="RMSE of the ensemble $(measure) for $(variable == "temp" ? "temperature" : "precipitation") on SSP119")
    hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$("ssp119").hdf5", "r")
    data = read(hfile, "rmse_$(measure)s_$(variable)_100")
    close(hfile)
    ext = (0., maximum(data))
    heatmap!(ax, lonvec2, latvec, data, colormap=:thermal, colorrange=ext)
    Colorbar(fig[1,2], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    display(fig)
    save("figs/ens_$(measure)_rmse_$(variable)_ssp119.png", fig)
end

#now test ks
ks = [x for x in 1:3]
variable = "pr" #temp/pr
for scenario in scenarios[2:end]
    # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

    hfile = h5open("data/temp_precip/vars_$(variable)_$(scenario)_50ens.hdf5", "r")
    true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
    close(hfile)

    wfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "cw") #let's make this also include means
    hfile = h5open("data/temp_precip/ens_vars_withpr_$(scenario).hdf5", "r") #this includes means
    for k in ks
        ens_means = variable == "temp" ? read(hfile, "ens_means_tas_k$(k)_new") : read(hfile, "ens_means_$(variable)_k$(k)_new")

        rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
        rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries

        write(wfile, "rmse_means_$(variable)_k$(k)_new", rmse_means)
        write(wfile, "rmse_means_time_$(variable)_k$(k)_new", rmse_means_time)
    end
    close(hfile)
    close(wfile)
end

measure = "mean" #mean only here
variable = "pr" #temp/pr
begin
    fig = Figure(resolution=(800,600))
    ax = Axis(fig[1,1], title="Average RMSE of the ensemble $(measure) for $(variable == "temp" ? "temperature" : "precipitation")", xlabel="Year", ylabel="RMSE")
    linestyles = [ :dot, :dash, :solid]
    for (j, scenario) in enumerate(scenarios[2:end])
        for (i, k) in enumerate(ks)
            hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "r")
            ser = read(hfile, "rmse_$(measure)s_time_$(variable)_k$(k)_new")
            lines!(ax,time_future, month_to_year_avg(ser), color=scenario_colors[scenario], alpha=0.6, linestyle=linestyles[i])
            close(hfile)
        end
    end

    elems = [LineElement(color=:black, linestyle=:dot),LineElement(color=:black, linestyle=:dash), LineElement(color=:black, linestyle=:solid)]
    elems_2 = [ LineElement(color=scenario_colors["ssp119"]), LineElement(color=scenario_colors["ssp245"]), LineElement(color=scenario_colors["ssp585"])]
    labels = ["k=1",  "k=2", "k=3", "SSP119", "SSP245", "SSP585"]
    axislegend(ax, [elems..., elems_2...], labels, position=:lt)
    display(fig)
    save("figs/ens_$(measure)_rmse_time_$(variable)_k_comparison.png", fig)
end

