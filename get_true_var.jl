using ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

M = 192
N = 96
# L1 = 1872 #for historical (generalize!)
# L2 = 1128 #for rcp (generalize!)
L1 = 1980
L2 = 1032

scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]
scenario = scenarios[parse(Int, ARGS[1])]
num_ens_members = 50
variables = ["tas", "pr"]
variable = variables[2]
folder_name = "temp_precip"

println("working on $(scenario)")

# all_gmts = zeros((num_ens_members, Int(L2/12))) #crucial! need to change this for rcp futures 
all_data = zeros((M,N, L2, num_ens_members))
for i in ProgressBar(1:num_ens_members)
    println("working on ensemble member $(i)")
    flush(stdout)
    file = "/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/$(variable)/r$(i)i1p1f1_$(scenario)_$(variable).nc"
    ts = ncData(file, variable) 
    # all_gmts[i,:] = get_gmt_list(ts)
    all_data[:,:,:,i] = ts.data[:,:,:]
end
true_var = var(all_data, dims=4)[:,:,:,1]
# true_ens_gmt = mean(all_gmts, dims=1)[:]
true_ens_mean = mean(all_data, dims=4)[:,:,:,1]

hfile = h5open("data/$(folder_name)/vars_$(variable)_$(scenario)_50ens.hdf5", "w")
write(hfile, "true_var", true_var)
write(hfile, "num_ens_members", num_ens_members)
# write(hfile, "true_ens_gmt", true_ens_gmt)
write(hfile, "true_ens_mean", true_ens_mean)
close(hfile)