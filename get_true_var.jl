using ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

M = 192
N = 96
L1 = 1872 #for historical (generalize!)
L2 = 1128 #for rcp (generalize!)

num_ens_members = 90

all_gmts = zeros((num_ens_members, Int(L2/12))) #crucial! need to change this for rcp futures 
all_data = zeros((M,N, L2, num_ens_members))
for i in 1:num_ens_members#ProgressBar(1:num_ens_members)
    println("working on ensemble member $(i)")
    flush(stdout)
    # file = "/net/fs06/d3/CMIP5/MPI-GE/RCP45/ts/ts_Amon_MPI-ESM_rcp45_r$(string(i, pad=3))i2005p3_200601-209912.nc"
    file = "/net/fs06/d3/CMIP5/MPI-GE/RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"
    ts = ncData(file, "ts") 
    # all_gmts[i,:] = get_gmt_list(ts)
    all_data[:,:,:,i] = ts.data[:,:,:]
end
# true_var = var(all_data, dims=4)[:,:,:,1]
# true_ens_gmt = mean(all_gmts, dims=1)[:]
true_ens_mean = mean(all_data, dims=4)

hfile = h5open("data/vars_rcp85_90ens.hdf5", "r+")
# write(hfile, "true_var", true_var)
# write(hfile, "num_ens_members", num_ens_members)
# write(hfile, "true_ens_gmt", true_ens_gmt)
write(hfile, "true_ens_mean", true_ens_mean)
close(hfile)