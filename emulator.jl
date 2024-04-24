# this is just an analysis script

using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list
file3 = "/net/fs06/d3/CMIP5/MPI-GE/RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(91, pad=3))i2005p3_200601-209912.nc" 
ts3 = ncData(file3, "ts") # actual data for comparison
gmt_list = get_gmt_list(ts3)
M, N, L = size(ts3.data)
latvec = ts3.latvec

M = 192
N = 96
L = 1128

# sim = emulate(gmt_list, mean_coefs, chol_coefs)
# sim = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
# newdata = back_to_data(sim, basis)
# simts = ncData(shape_data(newdata, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
# sim_gmt = get_gmt_list(simts)

#split the data
# newts = newdata[1:M*N,:]
# newtp = newdata[M*N+1:end, :]

#### correct way of comparing rmse -- only needs chol_coefs as the 'emulator' component! 

hfile = h5open("data/vars_rcp26_90ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
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
    L = 1128
    hfile = h5open("data/gaussian_emulator_rcp85_$(d)d.hdf5", "r")
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
# # ext_means = (0., maximum(rmse_means_20)) #this is FIXED based on teh rcp8.5 compariosn
# rmse_means_60 = sqrt.(sum((true_ens_mean.-ens_means_60).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])
# rmse_means_100 = sqrt.(sum((true_ens_mean.-ens_means_100).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])
# rmse_means_200 = sqrt.(sum((true_ens_mean.-ens_means_200).^2, dims=3)[:,:,1]./size(true_ens_mean)[3])

# rmse_means_time_20 = sqrt.(weighted_avg((true_ens_mean.-ens_means_20).^2, latvec)) #spatial average rmse (shaped as a timeseries)
# rmse_means_time_60 = sqrt.(weighted_avg((true_ens_mean.-ens_means_60).^2, latvec))
# rmse_means_time_100 = sqrt.(weighted_avg((true_ens_mean.-ens_means_100).^2, latvec))
# rmse_means_time_200 = sqrt.(weighted_avg((true_ens_mean.-ens_means_200).^2, latvec))

# rmse_mean_numbers = zeros((4, 3))
# all_rmse_means_time = zeros((4, 3, 1128))
# hfile = h5open("data/ens_means_rmse_comparison.hdf5", "w")
# write(hfile, "rmse_mean_numbers", rmse_mean_numbers)
# write(hfile, "all_rmse_means_time", all_rmse_means_time)
# close(hfile)
######################################################################################################################

# ens_vars_20 = get_ens_vars(20, true_ens_gmt)
# ens_vars_60 = get_ens_vars(60, true_ens_gmt)
# ens_vars_100 = get_ens_vars(100, true_ens_gmt)
# ens_vars_200 = get_ens_vars(200, true_ens_gmt)

# hfile = h5open("data/ens_vars_rcp26.hdf5", "r+")
# write(hfile, "ens_vars_20", ens_vars_20)
# write(hfile, "ens_vars_60", ens_vars_60)
# write(hfile, "ens_vars_100", ens_vars_100)
# write(hfile, "ens_vars_200", ens_vars_200)
# close(hfile)

hfile = h5open("data/ens_vars_rcp26.hdf5", "r")
ens_vars_20 = read(hfile, "ens_vars_20")
ens_vars_60 = read(hfile, "ens_vars_60")
ens_vars_100 = read(hfile, "ens_vars_100")
ens_vars_200 = read(hfile, "ens_vars_200")
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
# ext_var_rmse = (0., maximum(rmse_20)) - fix this based on the rcp8.5 comparison
rmse_60 = sqrt.(sum((true_var.-ens_vars_60).^2, dims=3)[:,:,1]./size(true_var)[3])
rmse_100 = sqrt.(sum((true_var.-ens_vars_100).^2, dims=3)[:,:,1]./size(true_var)[3])
rmse_200 = sqrt.(sum((true_var.-ens_vars_200).^2, dims=3)[:,:,1]./size(true_var)[3])

rmse_time_20 = sqrt.(weighted_avg((true_var.-ens_vars_20).^2, latvec)) #spatial average rmse (shaped as a timeseries
rmse_time_60 = sqrt.(weighted_avg((true_var.-ens_vars_60).^2, latvec))
rmse_time_100 = sqrt.(weighted_avg((true_var.-ens_vars_100).^2, latvec))
rmse_time_200 = sqrt.(weighted_avg((true_var.-ens_vars_200).^2, latvec))

# rcp26_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# rcp45_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)
# rcp85_time_avgs = (rmse_time_20, rmse_time_60, rmse_time_100, rmse_time_200)

# rmse_var_numbers = zeros((4, 3))
# all_rmse_vars_time = zeros((4, 3, 1128))
# hfile = h5open("data/ens_vars_rmse_comparison.hdf5", "w")
# write(hfile, "rmse_var_numbers", rmse_var_numbers)
# write(hfile, "all_rmse_vars_time", all_rmse_vars_time)
# close(hfile)

hfile = h5open("data/ens_means_rmse_comparison.hdf5", "r+")
rmse_mean_numbers = read(hfile, "rmse_mean_numbers")
close(hfile)


begin
    fig = Figure(resolution=(2000,1500))
    ax = Axis(fig[1,1], title="20 modes")
    # heatmap!(ax, rmse_20, colorrange=extrema(rmse))
    heatmap!(ax, rmse_means_20, colorrange=ext_means)
    ax = Axis(fig[1,2], title="60 modes")
    # heatmap!(ax, rmse_60,  colorrange=extrema(rmse))
    heatmap!(ax, rmse_means_60,  colorrange=ext_means)
    ax = Axis(fig[2,1], title="100 modes")
    # heatmap!(ax, rmse_100, colorrange=extrema(rmse))
    heatmap!(ax, rmse_means_100, colorrange=ext_means)
    ax = Axis(fig[2,2], title="200 modes")
    # heatmap!(ax, rmse_200, colorrange=extrema(rmse))
    heatmap!(ax, rmse_means_200, colorrange=ext_means)
    Colorbar(fig[1,3], colorrange=ext)
    display(fig)
    # save("figs/ens_mean_rmse_rcp26.png", fig)
end


begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1], title="average RMSE of the variance", xlabel="years since 2005", ylabel="RMSE")
    linestyles = [ :dot, :dashdot, :dash, :solid,]
    for (i, ser) in enumerate(rcp26_time_avgs)
        rmse_var_numbers[i, 1] = mean(ser)
        all_rmse_vars_time[i, 1, :] = ser
        lines!(ax, month_to_year_avg(ser), color=:orange, alpha=0.6, linestyle=linestyles[i])
    end
    for (i, ser) in enumerate(rcp45_time_avgs)
        rmse_var_numbers[i, 2] = mean(ser)
        all_rmse_vars_time[i, 2, :] = ser
        lines!(ax, month_to_year_avg(ser), color=:green, alpha=0.6, linestyle=linestyles[i])
    end
    for (i, ser) in enumerate(rcp85_time_avgs)
        rmse_var_numbers[i, 3] = mean(ser)
        all_rmse_vars_time[i, 3, :] = ser
        lines!(ax, month_to_year_avg(ser), color=:blue, alpha=0.6, linestyle=linestyles[i])
    end
    # lines!(ax, month_to_year_avg(rmse_time), label="20 modes")
    # lines!(ax, month_to_year_avg(rmse_time_60), label="60 modes")
    # lines!(ax, month_to_year_avg(rmse_time_100), label="100 modes")
    # lines!(ax, month_to_year_avg(rmse_time_200), label="200 modes")

    elems = [LineElement(color=:black, linestyle=:dot), LineElement(color=:black, linestyle=:dashdot), LineElement(color=:black, linestyle=:dash), LineElement(color=:black, linestyle=:solid)]
    elems_2 = [LineElement(color=:orange), LineElement(color=:green), LineElement(color=:blue)]
    labels = ["20 modes", "60 modes", "100 modes", "200 modes", "RCP2.6", "RCP4.5", "RCP8.5"]
    axislegend(ax, [elems..., elems_2...], labels, position=:lb)
    display(fig)
    # save("figs/ens_var_rmse_time_rcp_comparison.png", fig)
end






####### compare gmt representation w some ensemble members
ens_mem = 10
gmts = zeros((ens_mem, Int(L1/12)))
gmts_cov = zeros((ens_mem, Int(L1/12)))
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
        lines!(ax, gmts[i, :], color=:orange, alpha=0.3) #orange = without covariance accounted for
    end
    lines!(ax, avg, label="emulator ensemble avg", color=:red)
    lines!(ax, avg_cov, label="emulator ensemble avg with cov", color=:blue)
    lines!(ax,  gmt_list, label="model", color=:black)
    axislegend(ax, position=:lt)
    display(fig)
    # save("figs/sim_gmt_compare_cov_rcp85.png", fig)
end

#testing the ensemble variance - WRONG, but could be interesting to check simulated var and prescribed cov-derived
ens_mem = 90
sim_ens_projts = zeros((d, L1, ens_mem))
for i in ProgressBar(1:ens_mem)
    sim_ens_projts[:,:,i] = emulate(gmt_list, mean_coefs, chol_coefs; no_cov=true)
end
sim_ens_var = var(sim_ens_projts, dims=3)[:,:,1]









###### test other statistics? (older!!)

#Nino3.4 index
ensoslice = simts(-150, -90, -5, 5)
eM, eN, eL = size(ensoslice.data)
# ext = extrema(ensoslice.data[:,:,1])
# for i in 1:100
#     begin
#         fig = Figure(resolution=(1000, 200))
#         ax = Axis(fig[1,1])
#         heatmap!(ax, ensoslice.data[:,:,i], colorrange=ext)
#         display(fig)
#     end
# end
ets = zeros((20, Int(L1/12)))
for i in ProgressBar(1:20)
    s = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
    data = back_to_data(s, basis)
    sts = ncData(shape_data(data, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)(-150, -90, -5, 5)
    ets[i,:] = get_gmt_list(sts) #i think this works for anything --> yearly weighted avg timeseries
end
avg = dropdims(mean(ets, dims=1),dims=1)


obsslice = ts3(-150, -90, -5, 5)
# obs_indices = dropdims(mean(obsslice.data, dims=[1,2]), dims=(1,2))
# sim_indices = dropdims(mean(ensoslice.data, dims=[1,2]), dims=(1,2)) # could make this a non-naive average/ but also eq is fine

obs_indices = get_gmt_list(obsslice)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
        for i in 1:20
        lines!(ax, ets[i, :], color=:orange, alpha=0.3)
    end
    lines!(ax, avg, label="emulator", color=:red)
    lines!(ax,  obs_indices, label="model", color=:black)
    axislegend(ax)
    display(fig)
    # save("figs/enso34_ensemble.png", fig)
end

# (-120, -110, 34, 42) "nevada"
#preceip in boston/cape cod
file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(50, pad=3))i2005p3_200601-209912.nc"
tp3 = ncData(file, "pr")

pts = zeros((3, length(tp3.timevec)))
for i in ProgressBar(1:3)
    s = emulate(gmt_list, mean_coefs, corrs, var_coefs)
    data = back_to_data(s, basis)[M*N+1:end, :]
    sts = ncData(shape_data(data, M, N, true), tp3.lonvec, tp3.latvec, tp3.timevec)(-120, -102, 34, 45)
    # pts[i,:] = get_gmt_list(sts;yearly=false) #i think this works for anything --> yearly weighted avg timeseries
    pts[i,:] = weighted_avg(sts)
end
avg = dropdims(mean(pts, dims=1),dims=1)

obsslice = tp3(-120, -102, 34, 45)
# obs_indices = dropdims(mean(obsslice.data, dims=[1,2]), dims=(1,2))
# sim_indices = dropdims(mean(ensoslice.data, dims=[1,2]), dims=(1,2)) # could make this a non-naive average/ but also eq is fine
# obs_avgd = [mean(obs_indices[i:min(i+12-1, end)]) for i in 1:12:length(obs_indices)]
# obs_avgd = get_gmt_list(obsslice; yearly=false)
obs_avgd = weighted_avg(obsslice)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    for i in 1:3
        lines!(ax, pts[i, 12:48], label="emulator trajectory $i")#, color=:orange, alpha=1)
    end
    # lines!(ax, avg, label="emulator", color=:red)
    lines!(ax,  obs_avgd[12:48], label="model", color=:black, linewidth=3)
    axislegend(ax)
    display(fig)
    # save("figs/west_precip_3_trajectories_3yrs.png", fig)
end

using StatsBase

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    for i in 1:3
        lines!(ax, autocov(pts[i, :]), label="emulator trajectory $i")#, color=:orange, alpha=1)
    end
    # lines!(ax, avg, label="emulator", color=:red)
    lines!(ax,  autocov(obs_avgd), label="model", color=:black)
    axislegend(ax)
    display(fig)
    # save("figs/west_precip_autocorr_monthly.png", fig)
end