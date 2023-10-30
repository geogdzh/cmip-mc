using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra
include("utils.jl")
include("eof_util.jl")
# include("simulator.jl")

# using JLD
# all_eofs_cat = JLD.load("data/all_eofs_historical.jld", "all_eofs_cat")
# # all_eofs = [all_eofs_cat[:,:,i] for i in 1:100]

M = 192
N = 96

# U1 = all_eofs_cat[:,:,1] #U for first ensemble memeber

# eof1 = shape_data(U1[:,1], M, N)

# heatmap(eof1)
# extrema(eof1)

# basis = U1[:,1:5]

#sample snapshot:
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail
ts = ncData(file, "ts")
v = ts.data[:,:, 1]

timeseries = []
for i in 1:500
    v = ts.data[:,:,i]
    push!(timeseries,projection(reshape_data(v), basis))
end
timeseries = cat(timeseries..., dims=2)

lines(timeseries[2,:], timeseries[3,:])

##
X = reshape_data(ts.data)
U, S, V = svd(X)

basis = U[:,1:5]

lines(V[:,2], V[:,3], V[:,4])
lines(V[1:500,1])
lines(V[1:500,2])
lines(V[1:500,1],V[1:500,2],V[1:500,3])

#let's test smth out:
begin
    fig = Figure(resolution=(1000,1000))
    ax = Axis(fig[1,1])
    lines!(ax,timeseries[2,:], timeseries[3,:], label="timeseries projection",color=:red, alpha=0.3)
    lines!(ax,V[1:500,2].*S[2].+100,V[1:500,3].*S[3], label="V*S (+offset)", color=:blue, alpha=0.4)
    axislegend(ax)
    display(fig)
end
save( "figs/projection_sanity_check.png", fig)

# ok so let's see a map??

snap = shape_data(U[:,1].*S[1], M,N)

mode1_recovered = S[1] * reshape(U[:, 1], (192*96, 1)) * reshape(V[:, 1], (1, 1872))
mode2_recovered = S[2] * reshape(U[:, 2], (192*96, 1)) * reshape(V[:, 2], (1, 1872))
mode3_recovered = S[3] * reshape(U[:, 3], (192*96, 1)) * reshape(V[:, 3], (1, 1872))


# heatmap(shape_data(U[:,1].*S[1], M,N))
heatmap(shape_data(mode1_recovered[:, 1], M, N))
ext = extrema(mode1_recovered[:, 1].+ mode2_recovered[:, 1].+mode3_recovered[:, 1])

for i in 1:24
    fig = Figure(resolution=(1600,1000))
    ax = Axis(fig[1,1])
    summed = shape_data(mode1_recovered[:, i], M, N) + shape_data(mode2_recovered[:, i], M, N) + shape_data(mode3_recovered[:, i], M, N)
    heatmap!(ax, summed, colorrange=ext)
    display(fig)
end


tmp = back_to_data([1,1,1], U, S, V)
heatmap(shape_data(tmp[:,1], M, N))
tmp2 = mode1_recovered.+ mode2_recovered.+mode3_recovered
heatmap(shape_data(tmp2[:,1], M, N))

