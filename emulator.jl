using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, ProgressBars, GLM, Distributions, DataFrames, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

# get basis

#first ensemble member of historical run
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
ts = ncData(file_head*file_tail, "ts")
M, N, L1 = size(ts.data)
X = reshape_data(ts.data)

phfile = "/net/fs06/d3/CMIP5/MPI-GE/historical/precip/pr_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
pr = ncData(phfile, "pr")
Xp = reshape_data(pr.data)

#first ens member of 
file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
ts2 = ncData(file, "ts")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

phfile2 = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
pr2 = ncData(phfile2, "pr")
Xp2 = reshape_data(pr2.data)

# fullX = hcat(X, X2)
# U, S, V = svd(fullX)
# d = 100
# basis = U[:,1:d]

# hfile = h5open("data/temp_basis.hdf5", "w")
# write(hfile, "basis", basis)
# close(hfile)

fullX = hcat(vcat(X, Xp), vcat(X2, Xp2))
U, S, V = svd(fullX)
d = 100
basis = U[:,1:d]
# hfile = h5open("data/temp-precip_basis.hdf5", "w")
# write(hfile, "basis", basis)
# close(hfile)

################################### get the linear fits
num_ens_members = 3
projts = zeros((d, (L1+L2)*num_ens_members))
hist_mean_temp = zeros((1, (L1+L2)*num_ens_members))

# projts = project_timeseries(ts.data, basis)
# projts = hcat(projts, project_timeseries(ts2.data, basis))
# hist_mean_temp = vcat(weighted_avg(ts), weighted_avg(ts2))

#add more ensemble members
for i in ProgressBar(1:num_ens_members)
    files =  ["/net/fs06/d3/CMIP5/MPI-GE/historical/ts/ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc",
    "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
    files_pr = ["/net/fs06/d3/CMIP5/MPI-GE/historical/precip/pr_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc",
    "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
    for j in 1:2
        tmps = ncData(files[j], "ts")
        tmpp = ncData(files_pr[j], "pr")
        ind1, ind2 = (0,0)
        if j ==1
            ind1 = (i-1)*(L1+L2)+1
            ind2 = ind1-1+L1
        else
            ind1 = (i-1)*(L1+L2)+L1+1
            ind2 = ind1-1+L2
        end
        # projts[:, ind1:ind2] = project_timeseries(tmps.data, basis)
        data = vcat(reshape_data(tmps.data), reshape_data(tmpp.data))
        projts[:, ind1:ind2] = project_timeseries(data, basis; reshaped=true)
        hist_mean_temp[1, ind1:ind2] = weighted_avg(tmps)
    end
end

# using HDF5
# hfile = h5open("data/training_data_with26.hdf5", "w")
# write(hfile, "projts", projts)
# write(hfile, "hist_mean_temp", hist_mean_temp)
# write(hfile, "num_ens_members", num_ens_members)
# close(hfile)

hist_year_temp = [mean(hist_mean_temp[i:min(i+12-1, end)]) for i in 1:12:length(hist_mean_temp)]

# orig_covs =  zeros((2*d, 2*d, 12))  #debugging
corrs = zeros((2*d, 2*d, 12))
for i in 1:12 #i is prev month
    proj1 = projts[:,i:12:end]
    proj2 = projts[:,i+1:12:end]
    proj12 = [vcat(proj1[:,j], proj2[:,j]) for j in 1:size(proj2)[2]]
    proj12 = hcat(proj12...)
    corrs[:,:,i] = cor(proj12; dims=2)
    # orig_covs[:,:,i] = cov(proj12; dims=2)
end

###

#get mean coefs! #GENERALIZE
mean_coefs = zeros((12, d, 2))
for i in 1:12
    y = projts[:,i:12:end] 
    for j in 1:d
        mean_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=y[j,:])))
    end
end

#get var coefs! #GENERALIZE
var_coefs = zeros((12, d, 2))
for i in 1:12
    y = projts[:,i:12:end] 
    for j in 1:d
        b, m = mean_coefs[i, j, :]
        fits = [m*x+b for x in hist_year_temp]
        vars = (y[j,:].-fits).^2
        var_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=vars)))
    end
end


######################

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

### ok let's test it out for real

#get a sample gmt list
file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(50, pad=3))i2005p3_200601-209912.nc"
ts3 = ncData(file, "ts") # actual data for comparison
gmt_list = get_gmt_list(ts3)

sim = emulate(gmt_list, mean_coefs, corrs, var_coefs)
newdata = back_to_data(sim, basis)
# simts = ncData(shape_data(newdata, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
# sim_gmt = get_gmt_list(simts)

#split the data
newts = newdata[1:M*N,:]
newtp = newdata[M*N+1:end, :]

#compare gmt representation w 20 ensemble members
gmts = zeros((20, length(sim_gmt)))
for i in ProgressBar(1:20)
    s = emulate(gmt_list, mean_coefs, corrs, var_coefs)
    data = back_to_data(s, basis)[1:M*N, :]
    sts = ncData(shape_data(data, M, N, true), ts3.lonvec, ts3.latvec, ts3.timevec)
    gmts[i,:] = get_gmt_list(sts)
end
avg = dropdims(mean(gmts, dims=1),dims=1)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    # lines!(ax,  sim_gmt, label="emulator")
    for i in 1:20
        lines!(ax, gmts[i, :], color=:orange, alpha=0.3)
    end
    lines!(ax, avg, label="emulator ensemble avg", color=:red)
    lines!(ax,  gmt_list, label="model", color=:black)
    axislegend(ax)
    display(fig)
    save("figs/sim_gmt_20_members_with_precip.png", fig)
end


###### test other statistics?

#Nino3.4 index
ensoslice = simts(-150, -90, -5, 5)
eM, eN, eL = size(ensoslice.data)
ext = extrema(ensoslice.data[:,:,1])
for i in 1:100
    begin
        fig = Figure(resolution=(1000, 200))
        ax = Axis(fig[1,1])
        heatmap!(ax, ensoslice.data[:,:,i], colorrange=ext)
        display(fig)
    end
end


obsslice = ts3(-150, -90, -5, 5)
obs_indices = dropdims(mean(obsslice.data, dims=[1,2]), dims=(1,2))
sim_indices = dropdims(mean(ensoslice.data, dims=[1,2]), dims=(1,2)) # could make this a non-naive average/ but also eq is fine


begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    
    # lines!(ax,  sim_gmt, label="emulator")
    # for i in 1:20
    #     lines!(ax, gmts[i, :], color=:orange, alpha=0.3)
    # end
    lines!(ax, sim_indices[720:end], label="emulator", color=:red)
    lines!(ax,  obs_indices[720:end], label="model", color=:black)
    axislegend(ax)
    display(fig)
    # save("figs/enso34_3yrs.png", fig)
end


#preceip in boston/cape cod
file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/precip/pr_Amon_MPI-ESM_rcp26_r$(string(50, pad=3))i2005p3_200601-209912.nc"
tp3 = ncData(file, "pr")

pts = zeros((20, length(sim_gmt)))
for i in ProgressBar(1:20)
    s = emulate(gmt_list, mean_coefs, corrs, var_coefs)
    data = back_to_data(s, basis)[M*N+1:end, :]
    sts = ncData(shape_data(data, M, N, true), tp3.lonvec, tp3.latvec, tp3.timevec)(-72, -70, 41, 43)
    pts[i,:] = get_gmt_list(sts) #i think this works for anything --> yearly weighted avg timeseries
end
avg = dropdims(mean(pts, dims=1),dims=1)

obsslice = tp3(-72, -70, 41, 43)
obs_indices = dropdims(mean(obsslice.data, dims=[1,2]), dims=(1,2))
# sim_indices = dropdims(mean(ensoslice.data, dims=[1,2]), dims=(1,2)) # could make this a non-naive average/ but also eq is fine
# obs_avgd = [mean(obs_indices[i:min(i+12-1, end)]) for i in 1:12:length(obs_indices)]
obs_avgd = get_gmt_list(obsslice)

begin
    fig = Figure(resolution=(1000,800))
    ax = Axis(fig[1,1])
    
    # lines!(ax,  sim_gmt, label="emulator")
    for i in 1:20
        lines!(ax, pts[i, :], color=:orange, alpha=0.3)
    end
    lines!(ax, avg, label="emulator", color=:red)
    lines!(ax,  obs_avgd, label="model", color=:black)
    axislegend(ax)
    display(fig)
    # save("figs/enso34_3yrs.png", fig)
end

