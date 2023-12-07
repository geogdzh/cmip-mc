using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, Clustering, ColorSchemes
include("utils.jl")
include("simulator.jl")
include("eof_util.jl")

#### Step 1: define reference

#first ensemble member of historical run
file_head = "/net/fs06/d3/CMIP5/MPI-GE/historical/ts/"
file_tail = "ts_Amon_MPI-ESM_historical_r$(string(1, pad=3))i1850p3_185001-200512.nc"
file = file_head*file_tail
ts = ncData(file, "ts")

#get EOF basis
X = reshape_data(ts.data)
U, S, V = svd(X)
basis = U[:,1:5]

#partition state space based on the trajectory we got the eofs from
# trajectory = [V[:,i].*S[i] for i in 1:5]
# trajectory = transpose(hcat(trajectory...))

trajectories = Dict()
for j in 1:12
    traj = [V[j:12:end, i].*S[i] for i in 1:5]
    traj = transpose(hcat(traj...))
    trajectories[j] = traj
end

#trajectories is a dict of month isolated points

# k=5
# clusters = kmeans(trajectory, k)
# centers = clusters.centers

k=5
center_dict = Dict()
for j in 1:12
    clusters = kmeans(trajectories[j], k)
    center_dict[j] = clusters.centers
end

#figures to check cluster centers:
pal = colorschemes[:twelvebitrainbow]
begin
    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1], title="dims 2 vs 3")
    lines!(ax,V[1:500,2].*S[2].+100,V[1:500,3].*S[3])
    #scatter!(ax, centers[2,:], centers[3,:], color=:red, markersize=20)
    for j in 1:12
        centers = center_dict[j]
        scatter!(ax, centers[2,:], centers[3,:], color=pal[j], markersize=20)
    end
    display(fig)
    # save("figs/month_clusters_2_3.png", fig)
end

begin
    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1], title="dims 3 vs 4")
    lines!(ax,V[1:500,3].*S[3].+100,V[1:500,4].*S[4])#, label="V*S (+offset)", color=:blue, alpha=0.4)
    # scatter!(ax, centers[3,:], centers[4,:], color=:red, markersize=20)
    for j in 1:12
        centers = center_dict[j]
        scatter!(ax, centers[3,:], centers[4,:], color=pal[j], markersize=20)
    end
    # axislegend(ax)
    display(fig)
    # save("figs/month_clusters_3_4.png", fig)
end

begin
    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1], title="dims 1 vs 2")
    lines!(ax,V[1:500,1].*S[1].+100,V[1:500,2].*S[2])#, label="V*S (+offset)", color=:blue, alpha=0.4)
    # scatter!(ax, centers[1,:], centers[2,:], color=:red, markersize=20)
    for j in 1:12
        centers = center_dict[j]
        scatter!(ax, centers[1,:], centers[2,:], color=pal[j], markersize=20)
    end
    display(fig)
    # save("figs/month_clusters_1_2.png", fig)
end

#### Step 2: define partition function

function k_assignment(x, centers)
    distance(v1, v2) = sqrt(sum((v1 .- v2) .^ 2))
    return argmin([distance(x, col) for col in eachcol(centers)])
end

### Step 3: create Markov chain

function stepper(d::ncData, refs, assignment_func) #rght now it's just monthly
    #returns a dictionary, keyed on month, of markov chains for this dataset
    #refs should be of the format to fit the assignment function
    mc = []
    for t in 1:Int(length(d.timevec))
        month = Dates.month(d.timevec[t])
        proj = projection(d.data[:,:,t], basis)
        class = assignment_func(proj, refs[month]) + 5*(month-1)
        push!(mc, class)
    end
    return mc
end

mc = stepper(ts, center_dict, k_assignment)

using MarkovChainHammer: perron_frobenius

pf = perron_frobenius(mc)

heatmap(pf)

#try cutting pieces out of the pf:

monthly_pfs = Dict()
for j in 1:12
    monthly_pfs[j] = pf[5*(j-1)+1:5*j, 5*(j-1)+1:5*j]
end


