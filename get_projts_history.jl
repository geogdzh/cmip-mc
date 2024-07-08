using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

# L1, L2 = 1872, 1128 #for CMIP5
L1, L2 = 1980, 1032 #for CMIP6
using_precip = true 
non_dim = false  
use_metrics = true
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_precip
    parent_folder = "metrics"
elseif use_metrics && !using_precip
    parent_folder = "temp_metrics"
end

############## load a basis
hfile = using_precip ? h5open("data/$(parent_folder)/temp_precip_basis_1000d.hdf5", "r") : h5open("data/only_temp/temp_basis_1000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
if non_dim
    temp_factor = read(hfile, "temp_factor")
    pr_factor = read(hfile, "pr_factor")
end
close(hfile)
d = 200
basis = basis[:, 1:d]

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]

############### 
for scenario in scenarios[2:end]
    println("working on $(scenario)")
    flush(stdout)
    num_ens_members = 50 # number of model runs used to train the emulator
    projts = scenario == "historical" ? zeros((d, (L1), num_ens_members)) : zeros((d, (L2), num_ens_members))
    ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12)))

    file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
    # file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/"
    errors = []
    for i in 1:num_ens_members#ProgressBar(1:num_ens_members)
        try
            println("working on ensemble member $(i)")
            flush(stdout)
            files = file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"
            files_pr = file_head*"$(scenario)/pr/r$(i)i1p1f1_$(scenario)_pr.nc"
            tmps = ncData(files, "tas")
            prs = using_precip ? ncData(files_pr, "pr") : nothing
            if using_precip
                data = non_dim ? vcat(reshape_data(tmps.data) ./ temp_factor , reshape_data(prs.data) ./ temp_factor) : vcat(reshape_data(tmps.data), reshape_data(prs.data))
                projts[:, :, i] =  project_timeseries(data, basis, reshaped=true)
            else
                projts[:, :, i] =  project_timeseries(tmps.data, basis)
            end
            ens_gmt[i, :] = get_gmt_list(tmps) 
        catch
            println("missing values in ensemble member $(i)")
            push!(errors, i)
            flush(stdout)
        end
    end
    for i in reverse(errors)
    projts = projts[:,:,1:end .!= errors[:]] # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run # adn for some reason error when no errors?
    ens_gmt = ens_gmt[1:end .!= errors[:], :]
    end
    num_ens_members = size(ens_gmt)[1]

    hfile = using_precip ? h5open("data/$(parent_folder)/projts_withpr_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w") : h5open("data/temp_only/projts_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "projts", projts)
    write(hfile, "ens_gmt", ens_gmt)
    write(hfile, "num_ens_members", num_ens_members)
    close(hfile)
end