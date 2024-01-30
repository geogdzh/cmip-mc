using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, ProgressBars, GLM, Distributions, DataFrames#, ColorSchemes
include("utils.jl")
include("eof_util.jl")

#first ensemble member of historical run
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail
ts = ncData(file, "ts")

X = reshape_data(ts.data)
U, S, V = svd(X)
basis = U[:,1:5] #using basis from just 1 ensemble member - FIX TO ALSO INCLUDE RCP26

projts = project_timeseries(ts.data, basis)

proj_jan = projts[:,1:12:end]
proj_feb = projts[:,2:12:end]
janfeb = [vcat(proj_jan[:,i],proj_feb[:,i]) for i in 1:size(proj_jan)[2]]

for i in ProgressBar(2:20)
    file_tail = "ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc"
    file = file_head*file_tail
    ts = ncData(file, "ts")
    projts = project_timeseries(ts.data, basis)
    proj_jan = projts[:,1:12:end]
    proj_feb = projts[:,2:12:end]
    for j in 1:size(proj_jan)[2]
        push!(janfeb, vcat(proj_jan[:,i],proj_feb[:,i]))
    end
end

corr = cor(janfeb; dims=2)
covv = cov(janfeb; dims=2)


#check for gaussianity
janfeb_arr = hcat(janfeb...)
janfeb = hcat(janfeb...)
plot(janfeb_arr[1,:], janfeb_arr[10,:], alpha=0.5)


mvfit = fit(MvNormal, janfeb_arr)
Σ = mvfit.Σ

num_samples = 1000
samples = rand(mvfit, num_samples)

begin
    fig = Figure(resolution=(1000,1200))
    for i in 2:5
        ax = Axis(fig[i-1,1], title="feb 1 vs jan $i")
        plot!(ax, samples[6,:], samples[i,:], alpha=0.3, color=:red)
        plot!(ax, janfeb_arr[6,:], janfeb_arr[i,:], alpha=0.3, color=:black)
    end
    for i in 2:5
        ax = Axis(fig[i-1,2], title="feb 1 vs feb $i")
        plot!(ax, samples[6,:], samples[i+5,:], alpha=0.3, color=:red)
        plot!(ax, janfeb_arr[6,:], janfeb_arr[i+5,:], alpha=0.3, color=:black)
    end
    display(fig)
end


#####
hist_mean_temp = weighted_avg(ts)
hist_year_temp = [mean(hist_mean_temp[i:min(i+12-1, end)]) for i in 1:12:length(hist_mean_temp)]

# running_mean_temp = [sum(hist_mean_temp[i-12+1:i])/12 for i in 13:length(hist_mean_temp)] #starts in jan 1851, so 13:end on nay other file
# projts[:,13:end]
# running_mean_temp


mean_coefs = zeros((12, 5, 2)) #change from 5 to be general!!
for i in 1:12
    y = projts[:,i:12:end] #using just one ensemble member rn in projts
    for j in 1:5
        mean_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=y[j,:])))
    end
end

#check it for the first mode of january - looks good enough
b, m = mean_coefs[1, 1, :]
begin
    fig = Figure(resolution=(1000,1200))
    ax = Axis(fig[1,1])
    lines!(ax, hist_year_temp, [m*x+b for x in hist_year_temp])
    scatter!(ax, hist_year_temp, projts[:,1:12:end][1,:])
    display(fig)
end

#
vars_arr = zeros((12,5,156))
var_coefs = zeros((12, 5, 2))
for i in 1:12
    y = projts[:,i:12:end] #using just one ensemble member rn in projts
    for j in 1:5
        b, m = mean_coefs[i, j, :]
        fits = [m*x+b for x in hist_year_temp]
        vars = (y[j,:].-fits).^2
        vars_arr[i, j, :] = vars
        var_coefs[i, j, :]  = coef(lm(@formula(y ~ x), DataFrame(x=hist_year_temp, y=vars)))
    end
end

#check for jan
b, m = var_coefs[1, 1, :]
begin
    fig = Figure(resolution=(1000,1200))
    ax = Axis(fig[1,1])
    lines!(ax, hist_year_temp, [m*x+b for x in hist_year_temp])
    scatter!(ax, hist_year_temp, vars_arr[1,1,:])
    display(fig)
end

###
corr #the jan/feb 10x10 correlation matrix (but really only need to generate 12 of them once)

# size(corrs) = (10, 10, 12)




mens = get_means(289., mean_coefs)
covish = covv[6:end, 6:end]
dist = MvNormal(mens[:,1], covish)