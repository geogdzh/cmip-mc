using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

############## load a basis
# hfile = h5open("data/temp-precip_basis.hdf5", "r")
hfile = h5open("data/temp_basis.hdf5", "r") #this basis is calculated from just one ens member #and we're assuming that it's teh same for all rcp scenarios
basis = read(hfile, "basis")
close(hfile)
d = size(basis)[2]

#### load training data
hfile = h5open("data/training_data_with85_20d_90ens.hdf5", "r")
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
# ens_gmt = mean(ens_gmt, dims=1)
# ens_projts = projts
#OLD VERSION
# corrs = get_corrs(ens_projts) #updated
# mean_coefs = get_mean_coefs(ens_projts, ens_gmt) #updated
# var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=true) #updated

# NEW
mean_coefs = get_mean_coefs(ens_projts, ens_gmt) #updated
covs = gmt_cov(ens_projts, ens_gmt)
var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=true) #updated
corrs = gmt_cov(ens_projts, ens_gmt; corrs=true)
chol_coefs, chols = get_chol_coefs(covs, ens_gmt; return_chols=true)


# hfile = h5open("data/gaussian_emulator_rcp85_20d.hdf5", "w") #SHOULD add ensembel size?
# write(hfile, "mean_coefs", mean_coefs)
# write(hfile, "chol_coefs", chol_coefs)
# write(hfile, "basis", basis)
# close(hfile)

hfile = h5open("data/gaussian_emulator_rcp85_20d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
close(hfile)
d = size(basis)[2]

#################### ok let's test it out for real

#get a sample gmt list
# file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(50, pad=3))i2005p3_200601-209912.nc"
# file = "data/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(91, pad=3))i2005p3_200601-209912.nc" # 50 for what i was doing before
file = "data/ts/RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(91, pad=3))i2005p3_200601-209912.nc" # 50 for what i was doing before #FIX THE FILESYSTEM
ts3 = ncData(file, "ts") # actual data for comparison
gmt_list = get_gmt_list(ts3)
M, N, L1 = size(ts3.data)

# sim = emulate(gmt_list, mean_coefs, chol_coefs)
# sim = emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
# newdata = back_to_data(sim, basis)
# simts = ncData(shape_data(newdata, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
# sim_gmt = get_gmt_list(simts)

#split the data
# newts = newdata[1:M*N,:]
# newtp = newdata[M*N+1:end, :]

#### this feels much more correct
jan_1_cov = get_cov(true_ens_gmt[1], chol_coefs)[:,:,1]
tmp = sum([jan_1_cov[i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]) 
ens_vars = zeros(M, N, L1)
for m in ProgressBar(1:Int(L1/12))
    for n in 1:12
        co = get_cov(true_ens_gmt[m], chol_coefs)[:,:,n]
        ens_vars[:,:,(m-1)*12+n] = shape_data(sum([co[i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]), M, N)
    end
end
# heatmap(shape_data(tmp, M, N))
heatmap(ens_vars[:,:,13])

all_gmts = zeros((num_ens_members, Int(L1/12)))
all_data = zeros((M,N, L1, num_ens_members))
for i in ProgressBar(1:num_ens_members)
    file = "data/ts/RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"
    ts = ncData(file, "ts") 
    # all_gmts[i,:] = get_gmt_list(ts)
    all_data[:,:,:,i] = ts.data[:,:,:]
end
true_var = var(all_data, dims=4)[:,:,:,1]
ext = extrema(true_var)
true_ens_gmt = mean(all_gmts, dims=1)[:]

# hfile = h5open("data/vars_rcp85_90ens.hdf5", "w")
# write(hfile, "true_var", true_var)
# write(hfile, "num_ens_members", num_ens_members)
# write(hfile, "true_ens_gmt", true_ens_gmt)
# close(hfile)

begin
    monthtime = 1
    fig = Figure(resolution=(2000,800))
    ax = Axis(fig[1,1], title="true model variance")
    heatmap!(ax, true_var[:,:,monthtime], colorrange=ext)
    ax = Axis(fig[1,2], title="emulator variance")
    heatmap!(ax, ens_vars[:,:,monthtime], colorrange=ext)
    Colorbar(fig[1,3], colorrange=ext)
    display(fig)
    # save("figs/ens_var_compare_rcp85_month1.png", fig)
end

rmse = sqrt.(sum((true_var.-ens_vars).^2, dims=3)[:,:,1]./size(true_var)[3])
rmse_time = sqrt.(sum((true_var.-ens_vars).^2, dims=(1,2))[1,1,:]./(size(true_var)[1]*size(true_var)[2]))

lines(rmse_time[12:12:end])

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1], title="RMSE")
    heatmap!(ax, rmse)
    Colorbar(fig[1,2], colorrange=extrema(rmse))
    display(fig)
    save("figs/ens_var_rmse_rcp85_20modes.png", fig)
end

#compare gmt representation w some ensemble members
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

#testing the ensemble now

ens_mem = 90
sim_ens_projts = zeros((d, L1, ens_mem))
for i in ProgressBar(1:ens_mem)
    sim_ens_projts[:,:,i] = emulate(gmt_list, mean_coefs, chol_coefs; no_cov=true)
end
sim_ens_var = var(sim_ens_projts, dims=3)[:,:,1]

true_ens_projts = ens_projts[:,end+1-1128:end,:] #the future portion of it #FIX the hardcoding
true_ens_var = var(true_ens_projts, dims=3)[:,:,1] #method 1
data_true  = #isnt this just the original data lol? it is. (or this would be a compressed version of it)
    back_to_data(true_ens_projts[:,:,1], basis) #method 2
# heatmap(true_ens_var'-sim_ens_var',colormap=:balance)
timeind = 1
total = basis*sim_ens_var[:,timeind]
heatmap(shape_data(total, M, N))
ext = extrema(total)
total_true = basis*true_ens_var[:,timeind]
heatmap(shape_data(total_true, M, N), colorrange=ext)

begin
    fig = Figure(resolution=(2000,800))
    ax = Axis(fig[1,1], title="sim")
    heatmap!(ax, shape_data(total, M, N), colorrange=ext)
    # ax = Axis(fig[1,2], title="original model")
    # heatmap!(ax, shape_data(total_true, M, N), colorrange=ext)
    Colorbar(fig[1,3], colorrange=ext)
    display(fig)
    # save("figs/ens_var_compare_rcp85_month1.png", fig) # so it's representing SOMETHING well enough but it certainly isnt the ensemble variance
end




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


#####
# hfile = h5open("data/gaussian_emulator.hdf5", "w")
# write(hfile, "mean_coefs", mean_coefs)
# write(hfile, "var_coefs", var_coefs)
# write(hfile, "corrs", corrs)
# write(hfile, "basis", basis)
# close(hfile)


sim_var = var(simts.data[:,:,1:12:end], dims=3)
sample_var = var(ts3.data[:,:,1:12:end], dims=3)
ext = extrema(sample_var)
begin
    fig = Figure(resolution=(800,1000))
    ax = Axis(fig[1,1], title="sim")
    heatmap!(ax, sim_var[:,:,1], colorrange=ext)
    ax = Axis(fig[2,1], title="sample")
    heatmap!(ax, sample_var[:,:,1], colorrange=ext)
    Colorbar(fig[1,2], colorrange=ext)
    Colorbar(fig[2,2], colorrange=ext)
    display(fig)
    # save("figs/one_sample_var_jan_rcp85.png", fig)
end