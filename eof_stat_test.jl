using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
include("simulator.jl")
include("eof_util.jl")

file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail

d = ncData(file, "ts")
M = length(d.lonvec)
N = length(d.latvec)
L = length(d.timevec)


#ok so we know how to get a basis and work in it:
#want to then reconstruct... let's start with global mean temp?
X = reshape_data(d.data)
U, S, V = svd(X)

basisU = U[:, 1:5]
basisS = S[1:5]
basisV = V[:, 1:5]

###
newfile_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP85/ts/"
newfile_tail = "ts_Amon_MPI-ESM_rcp85_r$(string(1, pad=3))i2005p3_200601-209912.nc"
newfile = newfile_head*newfile_tail
newd = ncData(newfile, "ts")
newL = length(newd.timevec)

# project the new data onto original space
projts = project_timeseries(newd.data, basisU)

# visualize it:
lines(projts[2,1:500], projts[3,1:500]) #looks good but sufficiently different

###
# reconstruct "alternate universes" from each projected point
dims = (M*N, L, newL)
#reconstructed_proj = Array{Float64}(undef, dims)       #this gives an ot of memory error
#so instead go one by one
global_mean_temp_reduced = zeros((60))
for t in 1:60
    flatdata = back_to_data(projts[:,t], basisU, basisS, basisV)
    mean_temp = [mean(flatdata[:,i]) for i in 1:60] #could have been more than 60, doesn't have to be the same number
    global_mean_temp_reduced .+= mean_temp
end
global_mean_temp_reduced = global_mean_temp_reduced ./ 10
global_mean_temp_reduced = global_mean_temp_reduced ./ 6

lines(global_mean_temp_reduced)

#comparison:
global_mean_temp = Vector{Float32}()
for t in 1:60
   mean_temp = mean(newd.data[:,:,t]) 
   push!(global_mean_temp, mean_temp)
end

begin
    fig = Figure(resolution=(1600,1000))
    ax = Axis(fig[1,1])
    lines!(global_mean_temp_reduced, label="reduced")  #newd.timevec[1:60] ##timevec handling is wrong!!
    lines!( global_mean_temp, label="observed")
    axislegend()
    display(fig)
end