using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

############## load a basis
# hfile = h5open("data/temp-precip_basis.hdf5", "r")
hfile = h5open("data/temp_basis_20d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = size(basis)[2]

#### load training data
hfile = h5open("data/training_data_with26_20d.hdf5", "r")
ens_projts = read(hfile, "projts")
ens_gmt = read(hfile, "ens_gmt")
num_ens_members = read(hfile, "num_ens_members")
close(hfile)
ens_gmt = mean(ens_gmt, dims=1)

# hfile = h5open("data/training_data_with26.hdf5", "r")
# projts = read(hfile, "projts")
# hist_mean_temp = read(hfile, "hist_mean_temp")
# num_ens_members = read(hfile, "num_ens_members")
# close(hfile)

# d = size(projts)[1]
# ens_projts = zeros((d, Int(size(projts)[2]/num_ens_members), num_ens_members))
# for i in 1:num_ens_members
#     ens_projts[:,:,i] = projts[:,(i-1)*Int(size(projts)[2]/num_ens_members)+1:i*Int(size(projts)[2]/num_ens_members)]
# end

# hist_year_temp = month_to_year_avg(hist_mean_temp)
# ens_gmt = zeros((num_ens_members, Int(length(hist_year_temp)/num_ens_members)))
# for i in 1:num_ens_members
#     ens_gmt[i,:] = hist_year_temp[(i-1)*Int(length(hist_year_temp)/num_ens_members)+1:i*Int(length(hist_year_temp)/num_ens_members)]
# end
# ens_gmt = mean(ens_gmt, dims=1)
#now we're workign with ens_projts and ens_gmt

###
ens_gmt = mean(ens_gmt, dims=1)
ens_projts = projts
#OLD VERSION
# corrs = get_corrs(ens_projts) #updated
# mean_coefs = get_mean_coefs(ens_projts, ens_gmt) #updated
# var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=true) #updated

#NEW
mean_coefs = get_mean_coefs(ens_projts, ens_gmt) #updated
covs = gmt_cov(ens_projts, ens_gmt)
var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=true) #updated
# chol_coefs = get_chol_coefs(covs, ens_gmt)
corrs = gmt_cov(ens_projts, ens_gmt; corrs=true)
chol_coefs, chols = get_chol_coefs(covs, ens_gmt; return_chols=true)


####test
num_buckets = 10
bit = Int(size(ens_projts)[2]/num_buckets)
bitcorrs = zeros((2*d, 2*d, 12, num_buckets))
for i in 1:num_buckets
    bitcorrs[:,:,:,i] = get_corrs(ens_projts[:,(i-1)*bit+1:i*bit,:])
end
begin
    fig = Figure(resolution=(2000, 2000))
    month = 1
    ax = Axis(fig[1,1])
    heatmap!(ax, bitcorrs[d+1:end,1:d,month,5].-bitcorrs[d+1:end,1:d,month,1], colorrange=(-1,1),colormap=:balance)
    display(fig)
end
begin
    fig = Figure(resolution=(2000, 2000))
    # month = 1
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i,j], title = "mode $i vs mode $j of successive month")
            for month in 1:12
                lines!(ax, [x for x in 1:num_buckets], [bitcorrs[d+i,j,month,x] for x in 1:num_buckets])
            end
        end
    end
    # save("figs/corrs_over_time_10_buckets.png", fig)
    display(fig)
end


#mean fits
begin
    fig = Figure(resolution=(2000, 1000))
    # ax = Axis(fig[1,1])
    mon = 1
    # mode = 100
    for mode in 1:10
        # ax = Axis(fig[mode,1])
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[:], ens_projts[mode, mon:12:end, i], markersize=7, alpha=0.5) 
        end
        ens_mean = mean(ens_projts[mode, mon:12:end, :], dims=2)
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], [mean_coefs[mon, mode, 2].*x.+mean_coefs[mon, mode, 1] for x in ens_gmt[:]], color=:black, linewidth=4)
    end
    save("figs/jan_mean_fits_10_modes_rcp45.png", fig)
    display(fig)
end 

#var fits
begin
    fig = Figure(resolution=(2000, 1000))
    # ax = Axis(fig[1,1])
    mon = 1
    # mode = 100
    for mode in 1:10
        # ax = Axis(fig[mode,1])
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="std for mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[:], sqrt.(vars[mon,mode,i,:]), markersize=7, alpha=0.5) 
        end
        ens_mean = sqrt.(mean(dropdims(vars[mon:12:end, mode, :, :], dims=1), dims=1))
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], sqrt.([var_coefs[mon, mode, 2].*x.+var_coefs[mon, mode, 1] for x in ens_gmt[:]]), color=:black, linewidth=4)
    end
    save("figs/jan_var_fits_10_modes_rcp45.png", fig)
    display(fig)
end 

#cov fits
begin
    fig = Figure(resolution=(1400, 1400))
    mon = 1
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i,j])
            x = hcat(fill(1., length(ens_gmt)), ens_gmt[:])
            y = covs[d+i,j,mon,:] #CHANGE HERE for corrs/covs
            scatter!(ax, ens_gmt[:], y, color=:orange, alpha=0.5) 
            fits = x \ y
            # print(typeof(fit))
            lines!(ax, ens_gmt[:], fits[2].*ens_gmt[:] .+ fits[1], color=:black, linewidth=3)
            # lines!(ax, ens_gmt[:], [chol_coefs[month, i, j, 2].*x.+chol_coefs[month, i, j, 1] for x in ens_gmt[:]], color=:black, linewidth=3)
        end
    end
    save("figs/jan_cov_fits_rcp45.png", fig)
    display(fig)
end

#chol fits
begin
    fig = Figure(resolution=(1000, 1000))
    month = 1
    for i in 1:10
        for j in 1:10
            ax = Axis(fig[i,j])
            scatter!(ax, ens_gmt[:], chols[i,j,month,:], color=:orange, alpha=0.5) 
            lines!(ax, ens_gmt[:], [chol_coefs[month, i, j, 2].*x.+chol_coefs[month, i, j, 1] for x in ens_gmt[:]], color=:black, linewidth=3)
        end
    end
    display(fig)
end


######################
#testing
# covs = get_cov(289, corrs, var_coefs)
# means = get_means(289, mean_coefs)

# gmt_list = [290.0, 290.2, 290.3]
# sample = emulate(gmt_list, mean_coefs, corrs, var_coefs)
# newdata = back_to_data(sample, basis)

# ext = extrema(newdata[:,1])
# for i in 1:36
#     begin
#         fig = Figure(resolution=(1000, 800))
#         ax = Axis(fig[1,1])
#         heatmap!(ax, shape_data(newdata[:,i], M, N), colorrange=ext)
#         display(fig)
#     end
# end

#################### ok let's test it out for real

#get a sample gmt list
# file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(50, pad=3))i2005p3_200601-209912.nc"
file = "data/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(91, pad=3))i2005p3_200601-209912.nc" # 50 for what i was doing before
ts3 = ncData(file, "ts") # actual data for comparison
gmt_list = get_gmt_list(ts3)
M, N, L1 = size(ts3.data)

# sim = emulate(gmt_list, mean_coefs, chol_coefs)
sim = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
newdata = back_to_data(sim, basis)
simts = ncData(shape_data(newdata, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
sim_gmt = get_gmt_list(simts)

#split the data
# newts = newdata[1:M*N,:]
# newtp = newdata[M*N+1:end, :]

#compare gmt representation w some ensemble members
ens_mem = 10
gmts = zeros((ens_mem, length(sim_gmt)))
gmts_cov = zeros((ens_mem, length(sim_gmt)))
for i in ProgressBar(1:ens_mem)
    s = emulate(gmt_list, mean_coefs, chol_coefs)
    # s = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
    data = back_to_data(s, basis)[1:M*N, :]
    sts = ncData(shape_data(data, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
    gmts_cov[i,:] = get_gmt_list(sts)
end
avg = dropdims(mean(gmts, dims=1),dims=1)
avg_cov = dropdims(mean(gmts_cov, dims=1),dims=1)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    # lines!(ax,  sim_gmt, label="emulator")
    for i in 1:ens_mem
        lines!(ax, gmts[i, :], color=:orange, alpha=0.3)
        lines!(ax, gmts_cov[i, :], color=:blue, alpha=0.3)
    end
    lines!(ax, avg, label="emulator ensemble avg", color=:red)
    lines!(ax,  gmt_list, label="model", color=:black)
    axislegend(ax)
    display(fig)
    # save("figs/sim_gmt_20_members_with_precip.png", fig)
end


###### test other statistics?

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
ets = zeros((20, length(sim_gmt)))
for i in ProgressBar(1:20)
    s = emulate(gmt_list, mean_coefs, corrs, var_coefs)
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


#####
# hfile = h5open("data/gaussian_emulator.hdf5", "w")
# write(hfile, "mean_coefs", mean_coefs)
# write(hfile, "var_coefs", var_coefs)
# write(hfile, "corrs", corrs)
# write(hfile, "basis", basis)
# close(hfile)