using NCDatasets, ProgressBars
ENV["GKSwstype"] = "nul" 
using Statistics, SixelTerm, LinearAlgebra  # needed for the GR backend on headless servers
using CairoMakie

data_directory = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts"

file = "ts_Amon_MPI-ESM_historical_r" * string(1, pad = 3) * "i1850p3_185001-200512.nc"
filepath = joinpath(data_directory, file)
ds = Dataset(filepath)
starting_index = 50 * 12 + 1 
ending_index = 150 * 12
temperature = ds["ts"][:, :, starting_index:ending_index]
close(ds)
for i in ProgressBar(2:100)
    file = "ts_Amon_MPI-ESM_historical_r" * string(1, pad = 3) * "i1850p3_185001-200512.nc"
    filepath = joinpath(data_directory, file)
    ds = Dataset(filepath)
    temperature .+= ds["ts"][:, :, starting_index:ending_index]
    close(ds)
end
temperature .*= 1/100

##
mean_temp = mean(temperature, dims = (1,2))[:]


M = size(temperature, 1)
N = size(temperature, 2)
L = size(temperature, 3)
X = reshape(temperature, (M * N, L))
# ignores correct inner product
U, S, V = svd(X)

U100 = copy(U)
V100 = copy(V)

fig = Figure() 
ax11 = Axis(fig[1,1])
heatmap!(ax11, reshape(U[:, 1], (M, N)))
ax12 = Axis(fig[1,2])
heatmap!(ax12, reshape(U[:, 2], (M, N)))
ax21 = Axis(fig[2, 1])
heatmap!(ax21, reshape(U[:, 3], (M, N)))
ax22 = Axis(fig[2, 2])
heatmap!(ax22, reshape(U[:, 4], (M, N)))
display(fig)
##
# A ≈ U * Diagonal(S) * V'
R = reshape(reshape(U[:, 1], (M*N, 1)) * S[1] * reshape(V[:, 1], (1, L)), (M, N, L))
mean_R = mean(R, dims = (1,2))[:]
##
fig = Figure() 
ax11 = Axis(fig[1,1])
lines!(ax11, mean_R, color = (:red, 0.5))
lines!(ax11, mean_temp, color = (:blue, 0.5))
display(fig)