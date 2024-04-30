using ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]

# file3 = "/net/fs06/d3/mgeo/CMIP6/interim/ssp585/tas/r1i1p1f1_ssp585_tas.nc"
# ts3 = ncData(file3, "tas") 
# gmt_list = get_gmt_list(ts3)
# M, N, L = size(ts3.data)
# latvec = ts3.latvec


function get_ens_vars(d, true_ens_gmt; get_means=false) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = h5open("data/temp_precip/gaussian_emulator_withpr_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    close(hfile)
    ens_vars_tas = zeros(M, N, L)
    ens_vars_pr = zeros(M, N, L)
    for m in ProgressBar(1:Int(L/12))
        println("working on year $(m)")
        if get_means
            for n in 1:12
                data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1], basis)
                ens_vars_tas[:,:,(m-1)*12+n] = shape_data(data[1:M*N,:], M, N)
                ens_vars_pr[:,:,(m-1)*12+n] = shape_data(data[M*N+1:end,:], M, N)
            end
        else
            co = get_cov(true_ens_gmt[m], chol_coefs) 
            for n in 1:12
                data = sum([co[:,:,n][i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d])
                ens_vars_tas[:,:,(m-1)*12+n] = shape_data(data[1:M*N,:], M, N)
                ens_vars_pr[:,:,(m-1)*12+n] = shape_data(data[M*N+1:end,:], M, N)
            end
        end
    end
    println("done")
    return ens_vars_tas, ens_vars_pr
end


for scenario in scenarios[3:end] #CHANGE BACK
    println("working on $(scenario)")
    flush(stdout)
    #get the true gmt
    hfile = h5open("data/temp_only/vars_$(scenario)_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
    # true_var = read(hfile, "true_var")
    num_ens_members = read(hfile, "num_ens_members")
    true_ens_gmt = read(hfile, "true_ens_gmt")
    # true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
    close(hfile)

    ens_vars_tas_20, ens_vars_pr_20 = get_ens_vars(20, true_ens_gmt)
    ens_means_tas_20, ens_means_pr_20 = get_ens_vars(20, true_ens_gmt; get_means=true)
    println("working on 20")
    flush(stdout)
    hfile = h5open("data/temp_precip/ens_vars_withpr_$(scenario).hdf5", "w")
    write(hfile, "ens_vars_tas_20", ens_vars_tas_20)
    write(hfile, "ens_vars_pr_20", ens_vars_pr_20)
    write(hfile, "ens_means_tas_20", ens_means_tas_20)
    write(hfile, "ens_means_pr_20", ens_means_pr_20)
    close(hfile)
    for d in [100]#, 200] #CHANGE BACK
        println("working on $(d)")
        flush(stdout)
        ens_vars_tas, ens_vars_pr = get_ens_vars(d, true_ens_gmt)
        ens_means_tas, ens_means_pr = get_ens_vars(d, true_ens_gmt; get_means=true)
        hfile = h5open("data/temp_precip/ens_vars_withpr_$(scenario).hdf5", "r+")
        write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
        write(hfile, "ens_vars_pr_$(d)", ens_vars_pr)
        write(hfile, "ens_means_tas_$(d)", ens_means_tas)
        write(hfile, "ens_means_pr_$(d)", ens_means_pr)
        close(hfile)
    end

end