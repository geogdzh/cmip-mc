using ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]

# file3 = "/net/fs06/d3/mgeo/CMIP6/interim/ssp585/tas/r1i1p1f1_ssp585_tas.nc"
# ts3 = ncData(file3, "tas") 
# gmt_list = get_gmt_list(ts3)
# M, N, L = size(ts3.data)
# latvec = ts3.latvec

using_precip = true 
non_dim = false  
use_metrics = true
if using_precip
    parent_folder = "temp_precip"
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_precip
    parent_folder = "metrics"
elseif use_metrics && !using_precip
    parent_folder = "temp_metrics"
end


function get_ens_vars(d, true_ens_gmt; get_means=false, k=1) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = using_precip ? h5open("data/$(parent_folder)/gaussian_emulator_withpr_ssp585_$(d)d.hdf5", "r") : h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs_$(k)")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    if use_metrics
        metric = read(hfile, "metric") 
    end
    close(hfile)
    ens_vars_tas = zeros(M, N, L)
    ens_vars_pr = zeros(M, N, L)
    for m in ProgressBar(1:Int(L/12))
        println("working on year $(m)")
        flush(stdout)
        if get_means
            for n in 1:12
                if k==1 #MAKE THIS SMARTER lol
                    data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1], basis)
                elseif k==2
                    data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1] .+ mean_coefs[n,:,3].*true_ens_gmt[m].^2, basis)
                elseif k==3
                    data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1] .+ mean_coefs[n,:,3].*true_ens_gmt[m].^2 .+ mean_coefs[n,:,4].*true_ens_gmt[m].^3, basis)
                end
                ens_vars_tas[:,:,(m-1)*12+n] = use_metrics ? shape_data(data[1:M*N,:], M, N) ./ sqrt.(metric) : shape_data(data[1:M*N,:], M, N)
                if using_precip
                    ens_vars_pr[:,:,(m-1)*12+n] = use_metrics ? shape_data(data[M*N+1:end,:], M, N) ./ sqrt.(metric) : shape_data(data[M*N+1:end,:], M, N)
                end
            end
        else
            co = get_cov(true_ens_gmt[m], chol_coefs) 
            for n in 1:12
                data = sum([co[:,:,n][i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d])
                ens_vars_tas[:,:,(m-1)*12+n] = use_metrics ? shape_data(data[1:M*N,:], M, N) ./ sqrt.(metric) : shape_data(data[1:M*N,:], M, N)
                if using_precip
                    ens_vars_pr[:,:,(m-1)*12+n] = use_metrics ? shape_data(data[M*N+1:end,:], M, N) ./ sqrt.(metric) : shape_data(data[M*N+1:end,:], M, N)
                end
            end
        end
    end
    println("done")
    return using_precip ? (ens_vars_tas, ens_vars_pr) : ens_vars_tas
end

# test the differnt values of n
for scenario in scenarios[2:end] 
    println("working on $(scenario)")
    flush(stdout)
    #get the true gmt
    hfile = h5open("data/temp_only/vars_$(scenario)_50ens.hdf5", "r") 
    num_ens_members = read(hfile, "num_ens_members")
    true_ens_gmt = read(hfile, "true_ens_gmt")
    close(hfile)

    for d in  [20]#[20,100]#[200] #CHANGE BACK
        println("working on $(d)")
        flush(stdout)
        if using_precip
            ens_vars_tas, ens_vars_pr = get_ens_vars(d, true_ens_gmt)
            ens_means_tas, ens_means_pr = get_ens_vars(d, true_ens_gmt; get_means=true)
            hfile = h5open("data/$(parent_folder)/ens_vars_withpr_$(scenario).hdf5", "cw")
            write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
            write(hfile, "ens_vars_pr_$(d)", ens_vars_pr)
            write(hfile, "ens_means_tas_$(d)", ens_means_tas)
            write(hfile, "ens_means_pr_$(d)", ens_means_pr)
            close(hfile)
        else
            ens_vars_tas = get_ens_vars(d, true_ens_gmt)
            ens_means_tas = get_ens_vars(d, true_ens_gmt; get_means=true)
            hfile = h5open("data/$(parent_folder)/ens_vars_$(scenario).hdf5", "cw")
            write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
            write(hfile, "ens_means_tas_$(d)", ens_means_tas)
            close(hfile)
        end
    end
end

# # test the differnt values of k
# for scenario in scenarios[2:end] 
#     println("working on $(scenario)")
#     flush(stdout)
#     hfile = h5open("data/temp_only/vars_$(scenario)_50ens.hdf5", "r") 
#     num_ens_members = read(hfile, "num_ens_members")
#     true_ens_gmt = read(hfile, "true_ens_gmt")
#     close(hfile)

#     d=100  #this was the error
#     for k in 1:2
#         ens_means_tas, ens_means_pr = get_ens_vars(d, true_ens_gmt; get_means=true, k=k)
#         hfile = h5open("data/$(parent_folder)/ens_vars_withpr_$(scenario).hdf5", "cw")
#         write(hfile, "ens_means_tas_k$(k)", ens_means_tas)
#         write(hfile, "ens_means_pr_k$(k)", ens_means_pr)
#         close(hfile)
#     end
# end
