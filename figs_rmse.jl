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


## generate statistics
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

#now ks
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

## plots etc 

function plot_rmse(ax, variable, measure, numbers; testing_k=false)
    linestyles = testing_k ? [ :solid, :dash, :dashdot] : [ :dot, :solid]
    for (j, scenario) in enumerate(scenarios[2:end])
        for (i, number) in enumerate(numbers)
            hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$(scenario).hdf5", "r")
            if testing_k
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_k$(number)_new")
            else    
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_$(number)")
            end
            lines!(ax, time_future, sqrt.(month_to_year_avg(ser)), color=scenario_colors[scenario], alpha=0.6, linestyle=linestyles[i])
            close(hfile)
        end
    end

    elems = testing_k ? [LineElement(color=:black, linestyle=:solid),LineElement(color=:black, linestyle=:dash), LineElement(color=:black, linestyle=:dashdot)] : [LineElement(color=:black, linestyle=:solid), LineElement(color=:black, linestyle=:dot)]
    elems_2 = [ LineElement(color=scenario_colors["ssp119"]), LineElement(color=scenario_colors["ssp245"]), LineElement(color=scenario_colors["ssp585"])]
    labels = testing_k ? ["k=1",  "k=2", "k=3"] : [ "100 modes",  "20 modes"]
    labels_2 = [ "SSP119", "SSP245", "SSP585"]
    axislegend(ax, [elems..., elems_2...], [labels..., labels_2...], position=(measure=="var" && variable=="temp" ? :lb : :lt))
end


numbers = [20, 100]
ks = [x for x in 1:3]

variable = "temp" #temp/pr
begin #this version looks less like shit
    fig = Figure(resolution=(1500,1000)) #
    lims = Dict("temp" => (0.4, 0.8), "pr" => (0.0018, 0.0030))

    measure = "mean"
    ax = Axis(fig[1,1:4], title="a) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "precipitation") for varied # of modes", xlabel="Year", ylabel="RMSE")
    ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, numbers)
    ax = Axis(fig[1,5:8], title="b) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "precipitation") for varied degree of fit", xlabel="Year")
    ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, ks; testing_k=true)
    

    ax = GeoAxis(fig[2,1:5], title="d) RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "precipitation") on SSP119")
    hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$("ssp119").hdf5", "r")
    begin
        data = sqrt.(read(hfile, "rmse_$(measure)s_$(variable)_100"))
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data, colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,6], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end
    measure = "var"

    ax = Axis(fig[1,9:12], title="c) Average RMSE of the ensemble $(measure=="var" ? "std" : measure) \n for $(variable == "temp" ? "temperature" : "precipitation") for varied # of modes", xlabel="Year")
    plot_rmse(ax, variable, measure, numbers)
    
    ax = GeoAxis(fig[2,7:11], title="e) RMSE of the ensemble $(measure=="var" ? "std" : measure) \n for $(variable == "temp" ? "temperature" : "precipitation") on SSP119")
    hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$("ssp119").hdf5", "r")
    begin
        data = sqrt.(read(hfile, "rmse_$(measure)s_$(variable)_100"))
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data, colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,12], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end
    colsize!(fig.layout, 12, Relative(1/11))
    fig
    save("figs/rmse_joint_$(variable).png", fig)
end 

begin
    f = Figure()

    subgl_left = GridLayout()
    subgl_left[1, 1:3] = [Axis(f) for j in 1:3]

    subgl_right = GridLayout()
    # subgl_right[1, 1:4] = [Axis(f) for i in 1:4]

    f.layout[1, 1] = subgl_left
    f.layout[2, 1] = subgl_right
    ax4 = Axis(subgl_right[1,4])
    
    # heatmap!(subgl_right[1,1],  ratio[:,:,1])
    # Colorbar(subgl_right[1,2], colorrange=(0,1))
    # heatmap!(subgl_right[1,3], ratio[:,:,1])
    Colorbar(subgl_right[1,4], height = Relative(2/4))
    hidedecorations!(ax4)
    f
end



####### gridlayout
# begin
#     fig = Figure() #resolution=(1500,1000)
#     top = GridLayout()
#     bottom = GridLayout()
#     # top[1, 1:3] = [Axis(fig) for j in 1:3]
#     f.layout[1, 1] = top
#     f.layout[2, 1] = bottom

#     measure = "mean"
#     top[1,1] = Axis(fig, title="Average RMSE of the ensemble $(measure) for $(variable == "temp" ? "temperature" : "precipitation")", xlabel="Year", ylabel="RMSE")
#     plot_rmse(top[1,1], variable, measure, numbers)
#     top[1,2] = Axis(fig)
#     plot_rmse(top[1,2], variable, measure, ks; testing_k=true)

#     bottom[1,1] = GeoAxis(fig, title="RMSE of the ensemble $(measure) for $(variable == "temp" ? "temperature" : "precipitation") on SSP119")
#     hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$("ssp119").hdf5", "r")
#     begin
#         data = sqrt.(read(hfile, "rmse_$(measure)s_$(variable)_100"))
#         close(hfile)
#         ext = (0., maximum(data))
#         heatmap!(bottom[1,1], lonvec2, latvec, data, colormap=:thermal, colorrange=ext)
#         Colorbar(bottom[1,2], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
#     end
#     measure = "var"

#     top[1,3] = Axis(fig, title="Average RMSE of the ensemble $(measure=="var" ? "std" : measure) for $(variable == "temp" ? "temperature" : "precipitation")", xlabel="Year", ylabel="RMSE")
#     plot_rmse(top[1,3], variable, measure, numbers)
    
#     bottom[1,3] = GeoAxis(fig, title="RMSE of the ensemble $(measure=="var" ? "std" : measure) for $(variable == "temp" ? "temperature" : "precipitation") on SSP119")
#     hfile = h5open("data/temp_precip/ens_vars_withpr_rmse_$("ssp119").hdf5", "r")
#     begin
#         data = sqrt.(read(hfile, "rmse_$(measure)s_$(variable)_100"))
#         close(hfile)
#         ext = (0., maximum(data))
#         heatmap!(bottom[1,3], lonvec2, latvec, data, colormap=:thermal, colorrange=ext)
#         Colorbar(bottom[1,4], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
#     end

#     fig
#     # save("figs/rmse_joint_$(variable).png", fig)
# end 

hfile = h5open("data/temp_precip/vars_$("pr")_$("ssp585")_50ens.hdf5", "r")
future_mean = read(hfile, "true_ens_mean")
future_std = sqrt.(read(hfile, "true_var"))
close(hfile)
hfile = h5open("data/temp_precip/vars_$("pr")_$("historical")_49ens.hdf5", "r")
past_mean = read(hfile, "true_ens_mean")
past_std = sqrt.(read(hfile, "true_var"))
close(hfile)

diff_ratio = (mean(future_mean, dims=3) .- mean(past_mean, dims=3)) ./ mean(past_std, dims=3)
heatmap(diff_ratio[:,:,1], colorrange=(0,2), colormap=:balance)
extrema(diff_ratio)

std_ratio = mean(future_std, dims=3) ./ mean(past_std, dims=3)
heatmap(std_ratio[:,:,1], colorrange=(0,2), colormap=:balance)