using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, ProgressBars, GLM, Distributions, DataFrames#, ColorSchemes
include("utils.jl")
include("eof_util.jl")


# get basis

#first ensemble member of historical run
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
ts = ncData(file_head*file_tail, "ts")
X = reshape_data(ts.data)

file = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
ts2 = ncData(file, "ts")
X2 = reshape_data(ts2.data)
M, N, L = size(ts.data)

fullX = hcat(X, X2)
U, S, V = svd(fullX)
d = 5
basis = U[:,1:d]

################################### get the linear fits

projts = project_timeseries(ts.data, basis)
projts = hcat(projts, project_timeseries(ts2.data, basis))

orig_covs =  zeros((2*d, 2*d, 12)) 
corrs = zeros((2*d, 2*d, 12))
for i in 1:12 #i is prev month
    proj1 = projts[:,i:12:end]
    proj2 = projts[:,i+1:12:end]
    proj12 = [vcat(proj1[:,j], proj2[:,j]) for j in 1:size(proj2)[2]]
    proj12 = hcat(proj12...)
    # could then augment it with more ensemble members
    corrs[:,:,i] = cor(proj12; dims=2)
    orig_covs[:,:,i] = cov(proj12; dims=2)
end

###
hist_mean_temp = vcat(weighted_avg(ts), weighted_avg(ts2))
hist_year_temp = [mean(hist_mean_temp[i:min(i+12-1, end)]) for i in 1:12:length(hist_mean_temp)]

#get mean coefs! #GENERALIZE
mean_coefs = zeros((12, d, 2))
for i in 1:12
    y = projts[:,i:12:end] #using just one ensemble member rn in projts
    for j in 1:5
        mean_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=y[j,:])))
    end
end

#get var coefs! #GENERALIZE
var_coefs = zeros((12, d, 2))
for i in 1:12
    y = projts[:,i:12:end] #using just one ensemble member rn in projts
    for j in 1:5
        b, m = mean_coefs[i, j, :]
        fits = [m*x+b for x in hist_year_temp]
        vars = (y[j,:].-fits).^2
        var_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=vars)))
    end
end


######################
include("emulator_util.jl")

covs = get_cov(289, corrs, var_coefs)
means = get_means(289, mean_coefs)

gmt_list = [290.0, 290.2, 290.3]
sample = emulate(gmt_list, mean_coefs, corrs, var_coefs)
newdata = back_to_data(sample, basis)

ext = extrema(newdata[:,1])
for i in 1:36
    begin
        fig = Figure(resolution=(1000, 800))
        ax = Axis(fig[1,1])
        heatmap!(ax, shape_data(newdata[:,i], M, N), colorrange=ext)
        display(fig)
    end
end

### ok let's test it out for really