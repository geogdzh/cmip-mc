using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

L1, L2 = 1980, 1032 #for CMIP6
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

############## load a basis
hfile = using_precip ? h5open("data/$(parent_folder)/temp_precip_basis_2000d.hdf5", "r") : h5open("data/$(parent_folder)/temp_basis_2000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
if non_dim
    temp_factor = read(hfile, "temp_factor")
    pr_factor = read(hfile, "pr_factor")
end
if use_metrics
    metric = read(hfile, "metric")
end
close(hfile)
d = parse(Int, ARGS[1])
basis = basis[:, 1:d]

# scenario = "ssp119"
scenario = "ssp585"

############### get training data for the linear fits 
num_ens_members = 50 # number of model runs used to train the emulator
projts = zeros((d, (L1+L2), num_ens_members)) 
ens_gmt = zeros((num_ens_members, Int((L1+L2)/12))) #GMT in any case

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
# file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/"
errors = []
for i in 1:num_ens_members#ProgressBar(1:num_ens_members)
    try
        println("working on ensemble member $(i)")
        flush(stdout)
        files = [ file_head*"historical/tas/r$(i)i1p1f1_historical_tas.nc",
            file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"]
        files_pr = [ file_head*"historical/pr/r$(i)i1p1f1_historical_pr.nc",
            file_head*"$(scenario)/pr/r$(i)i1p1f1_$(scenario)_pr.nc"]
        tmps = ncData(files[1], "tas")
        prs = using_precip ? ncData(files_pr[1], "pr") : nothing
        if using_precip
            data1 = tmps.data
            data2 = prs.data
            if non_dim
                data1 = data1 ./ temp_factor
                data2 = data2 ./ pr_factor
            end
            # if use_metrics
            #     data1 = sqrt.(metric) .* data1
            #     data2 = sqrt.(metric) .* data2
            # end
            data = vcat(reshape_data(data1), reshape_data(data2))
            projts[:, 1:L1, i] =  project_timeseries(data, basis, reshaped=true)
        else
            # projts[:, 1:L1, i] =  use_metrics ? project_timeseries((tmps.data .* sqrt.(metric)), basis) : project_timeseries(tmps.data, basis)
            projts[:, 1:L1, i] = project_timeseries(tmps.data, basis)
        end
        ens_gmt[i, 1:Int(L1/12)] = get_gmt_list(tmps) # i had a use metrics here before but pretty sure that's wrong
        tmps = ncData(files[2], "tas")
        prs = using_precip ? ncData(files_pr[2], "pr") : nothing
        if using_precip
            data1 = tmps.data
            data2 = prs.data
            if non_dim
                data1 = data1 ./ temp_factor
                data2 = data2 ./ pr_factor
            end
            # if use_metrics
            #     data1 = sqrt.(metric) .* data1
            #     data2 = sqrt.(metric) .* data2
            # end
            data = vcat(reshape_data(data1), reshape_data(data2))
            projts[:, L1+1:end, i] = project_timeseries(data, basis, reshaped=true)
        else
            # projts[:, L1+1:end, i] = use_metrics ? project_timeseries((tmps.data .* sqrt.(metric)), basis) : project_timeseries(tmps.data, basis)
            projts[:, L1+1:end, i] =  project_timeseries(tmps.data, basis)
        end
        ens_gmt[i, Int(L1/12)+1:end] = get_gmt_list(tmps) # i had a use metrics here before but pretty sure that's wrong
    catch
        println("missing values in ensemble member $(i)")
        push!(errors, i)
        flush(stdout)
    end
end
# for i in reverse(errors)
projts = projts[:,:,1:end .!= errors[:]] # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run
ens_gmt = ens_gmt[1:end .!= errors[:], :]
# end
num_ens_members = size(ens_gmt)[1]

hfile = using_precip ? h5open("data/$(parent_folder)/training_data_withpr_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w") : h5open("data/$(parent_folder)/training_data_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
write(hfile, "projts", projts)
write(hfile, "ens_gmt", ens_gmt)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)

############## dev

# hfile = h5open("data/$(parent_folder)/training_data_$(scenario)_$(d)d_49ens.hdf5", "r")
# projts = read(hfile, "projts")
# ens_gmt = read(hfile, "ens_gmt")
# close(hfile)

# begin
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     for i in 1:49
#         lines!(ax, projts[1,1:12:end,i], color=:blue, label="tas")
#     end
#     display(fig)
# end
