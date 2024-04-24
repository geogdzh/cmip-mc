using ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]

file3 = "/net/fs06/d3/mgeo/CMIP6/interim/ssp585/tas/r$(1)i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas") 
gmt_list = get_gmt_list(ts3)
M, N, L = size(ts3.data)
latvec = ts3.latvec


function get_ens_vars(d, true_ens_gmt; get_means=false) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = h5open("data/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    close(hfile)
    ens_vars = zeros(M, N, L)
    for m in ProgressBar(1:Int(L/12))
        for n in 1:12
            if get_means
                ens_vars[:,:,(m-1)*12+n] = shape_data(back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1], basis), M, N)
            else
                co = get_cov(true_ens_gmt[m], chol_coefs)[:,:,n] #there is maybe a more efficient way to do this?
                ens_vars[:,:,(m-1)*12+n] = shape_data(sum([co[i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]), M, N)
            end
        end
    end
    return ens_vars
end


for scenario in scenarios[2:end]
    println("working on $(scenario)")
    flush(stdout)
    #point of comparison:
    hfile = h5open("data/vars_$(scenario)_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
    # true_var = read(hfile, "true_var")
    num_ens_members = read(hfile, "num_ens_members")
    true_ens_gmt = read(hfile, "true_ens_gmt")
    # true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
    close(hfile)

    ens_vars_20 = get_ens_vars(20, true_ens_gmt)
    println("working on 20")
    flush(stdout)
    hfile = h5open("data/ens_vars_$(scenario).hdf5", "w")
    write(hfile, "ens_vars_20", ens_vars_20)
    close(hfile)
    for d in [60, 100, 200]
        println("working on $(d)")
        flush(stdout)
        ens_vars = get_ens_vars(d, true_ens_gmt)
        hfile = h5open("data/ens_vars_$(scenario).hdf5", "r+")
        write(hfile, "ens_vars_$(d)d", ens_vars)
        close(hfile)
    end

end