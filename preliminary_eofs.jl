ENV["GKSwstype"] = "nul" 
using Statistics, SixelTerm, LinearAlgebra  # needed for the GR backend on headless servers
using CairoMakie
scatter(randn(10))

using NCDatasets
data_directory = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts"
file = "ts_Amon_MPI-ESM_historical_r" * string(100) * "i1850p3_185001-200512.nc"

filepath = joinpath(data_directory, file)
ds = Dataset(filepath)

time = ds["time"] 
# starting year 1850, ending year 2005 
# time[601] this for 20 years

starting_index = 50 * 12 + 1 
ending_index = starting_index + 20 * 12
temperature = ds["ts"][:, :, starting_index:ending_index]

m100 = mean(temperature)
close(ds)
heatmap(temperature[:, :, 1])

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
data_directory = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts"
file = "ts_Amon_MPI-ESM_historical_r" * "001" * "i1850p3_185001-200512.nc"
filepath = joinpath(data_directory, file)
ds = Dataset(filepath)

time = ds["time"] 
# starting year 1850, ending year 2005 
# time[601] this for 20 years

# starting_index = 601 
# ending_index = 601 + 20 * 12
temperature = ds["ts"][:, :, starting_index:ending_index]

close(ds)
heatmap(temperature[:, :, 1])
m001 = mean(temperature)


M = size(temperature, 1)
N = size(temperature, 2)
L = size(temperature, 3)
X = reshape(temperature, (M * N, L))
# ignores correct inner product
U, S, V = svd(X)


U001 = copy(U)
V001 = copy(V)


fig = Figure(resolution = (1200, 600)) 
ax11 = Axis(fig[1,1])
heatmap!(ax11, reshape(U100[:, 1], (M, N)), colorrange = extrema(U001[:, 1]))
ax12 = Axis(fig[1,2])
heatmap!(ax12, reshape(U100[:, 2], (M, N)), colorrange = extrema(U001[:, 2]))
ax21 = Axis(fig[2, 1])
heatmap!(ax21, reshape(U100[:, 3], (M, N)), colorrange = extrema(U001[:, 3]))
ax22 = Axis(fig[2, 2])
heatmap!(ax22, reshape(U100[:, 4], (M, N)), colorrange = extrema(U001[:, 4]))

ax13 = Axis(fig[1,3])
heatmap!(ax13, reshape(U001[:, 1], (M, N)), colorrange = extrema(U001[:, 1]))
ax14 = Axis(fig[1,4])
heatmap!(ax14, reshape(U001[:, 2], (M, N)), colorrange = extrema(U001[:, 2]))
ax23 = Axis(fig[2, 3])
heatmap!(ax23, reshape(U001[:, 3], (M, N)), colorrange = extrema(U001[:, 3]))
ax24 = Axis(fig[2, 4])
heatmap!(ax24, reshape(U001[:, 4], (M, N)), colorrange = extrema(U001[:, 4]))

display(fig)
save("eofs.png", fig, px_per_unit = 2)

##

fig = Figure(resolution = (1200, 600)) 
mode_shift = 4
ax11 = Axis(fig[1,1])
heatmap!(ax11, reshape(U100[:, 1+ mode_shift], (M, N)), colorrange = extrema(U001[:, 1+ mode_shift]))
ax12 = Axis(fig[1,2])
heatmap!(ax12, reshape(U100[:, 2+ mode_shift], (M, N)), colorrange = extrema(U001[:, 2+ mode_shift]))
ax21 = Axis(fig[2, 1])
heatmap!(ax21, reshape(U100[:, 3+ mode_shift], (M, N)), colorrange = extrema(U001[:, 3+ mode_shift]))
ax22 = Axis(fig[2, 2])
heatmap!(ax22, reshape(U100[:, 4+ mode_shift], (M, N)), colorrange = extrema(U001[:, 4+ mode_shift]))

ax13 = Axis(fig[1,3])
heatmap!(ax13, reshape(-U001[:, 1+ mode_shift], (M, N)), colorrange = extrema(U001[:, 1+ mode_shift]))
ax14 = Axis(fig[1,4])
heatmap!(ax14, reshape(-U001[:, 2+ mode_shift], (M, N)), colorrange = extrema(U001[:, 2+ mode_shift]))
ax23 = Axis(fig[2, 3])
heatmap!(ax23, reshape(-U001[:, 3+ mode_shift], (M, N)), colorrange = extrema(U001[:, 3+ mode_shift]))
ax24 = Axis(fig[2, 4])
heatmap!(ax24, reshape(U001[:, 4+ mode_shift], (M, N)), colorrange = extrema(U001[:, 4+ mode_shift]))

display(fig)
save("higher_modes_eofs.png", fig, px_per_unit = 2)

##

diffs = [norm(abs.(U100[:, i]) - abs.(U001[:, i])) for i in 1:100]

fig = Figure(resolution = (600, 600)) 
ax11 = Axis(fig[1,1])
lm1 = maximum(U100[:, 1])
heatmap!(ax11, reshape(U100[:, 1] - U001[:, 1], (M, N)), colormap = :balance, colorrange = (-lm1, lm1))
ax12 = Axis(fig[1,2])
lm1 = maximum(U100[:, 2])
heatmap!(ax12, reshape(U100[:, 2]- U001[:, 2], (M, N)), colormap = :balance, colorrange = (-lm1, lm1))
ax21 = Axis(fig[2, 1])
lm1 = maximum(U100[:, 3])
heatmap!(ax21, reshape(U100[:, 3]- U001[:, 3], (M, N)), colormap = :balance, colorrange = (-lm1, lm1))
ax22 = Axis(fig[2, 2])
lm1 = maximum(U100[:, 4])
heatmap!(ax22, reshape(U100[:, 4]- U001[:, 4], (M, N)), colormap = :balance, colorrange = (-lm1, lm1))
display(fig)
save("eofs_diffs.png", fig, px_per_unit = 2)