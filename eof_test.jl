using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
include("simulator.jl")

file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail

# file_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP26/ts/"
# file_tail = "ts_Amon_MPI-ESM_rcp26_r$(string(1, pad=3))i2005p3_200601-209912.nc"
# file = file_head*file_tail

# file_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP45/ts/"
# file_tail = "ts_Amon_MPI-ESM_rcp45_r$(string(1, pad=3))i2005p3_200601-209912.nc"
# file = file_head*file_tail

# file_head = "/net/fs06/d3/CMIP5/MPI-GE/RCP85/ts/"
# file_tail = "ts_Amon_MPI-ESM_rcp85_r$(string(1, pad=3))i2005p3_200601-209912.nc"
# file = file_head*file_tail

d = ncData(file, "ts")
M = length(d.lonvec)
N = length(d.latvec)


# first: let's check the EOFs of all the different ensemble members
# get all EOFs, visualize them as the variance across all ensemble members
#(will then do the same for all the scenarios / across scenarios)

all_eofs = []

for i in 1:100
    # file_tail = "ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc"
    file_tail = "ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"
    file = file_head*file_tail
    d = ncData(file, "ts")
    U = get_eof(d)
    push!(all_eofs, U)
end

using JLD
all_eofs_cat = cat(all_eofs...,dims=3)
# JLD.save("all_eofs_historical.jld", "all_eofs_cat", all_eofs_cat)
# JLD.save("data/all_eofs_rcp26.jld", "all_eofs_cat", all_eofs_cat)
# JLD.save("data/all_eofs_rcp45.jld", "all_eofs_cat", all_eofs_cat)
# JLD.save("data/all_eofs_rcp85.jld", "all_eofs_cat", all_eofs_cat)

all_eofs_cat = JLD.load("data/all_eofs_historical.jld", "all_eofs_cat")
all_eofs = [all_eofs_cat[:,:,i] for i in 1:100]


# check xth eof:
x = 5

begin
    eof1s = []
    for U in all_eofs
        eof1 = U[:, x]
        map = shape_data(eof1, M, N)
        push!(eof1s, map)
    end
    eof1_ag = cat(eof1s..., dims=3)
    eof1_var = var(eof1_ag, dims=3)[:, :, 1]
end
begin
    fig = Figure(resolution = (1200, 600))
    ax = Axis(fig[1,1]) 
    heatmap!(ax,eof1_var)
    Colorbar(fig[1,2])
    display(fig)

end

heatmap(shape_data(all_eofs[4][:,2], M, N))
maximum(all_eofs[1])

#compare first five eof's and their variances
fig = Figure(resolution=(6000,1200))
lims = extrema(all_eofs[1][:,5]) #(-0.02758066f0, 0.031166947f0)
std_lims = extrema(sqrt.(eof1_var)) #(0.00011589288f0, 0.03073324f0)
for x in 1:5

    ax1 = Axis(fig[1,x])
    heatmap!(ax1, shape_data(all_eofs[1][:,x], M, N), colorrange=lims)

    eof1s = []
    for U in all_eofs
        eof1 = U[:, x]
        map = shape_data(eof1, M, N)
        push!(eof1s, map)
    end
    eof1_ag = cat(eof1s..., dims=3)
    eof1_std = std(eof1_ag, dims=3)[:, :, 1]

    ax2 = Axis(fig[2,x])
    heatmap!(ax2, eof1_std, colorrange=std_lims)

end
display(fig)
save("./figs/eof_std_historical.png", fig, px_per_unit=2)


#now, the other question is to compare across scenarios - 
#maybe a map of avg difference across 