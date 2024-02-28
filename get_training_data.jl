using HDF5
include("utils.jl")
include("eof_util.jl")

L1 = 1872 #for historical (generalize!)
L2 = 1128 #for rcp26 (generalize!)

############## load a basis
# hfile = h5open("data/temp-precip_basis.hdf5", "r")
hfile = h5open("data/temp_basis_20d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = size(basis)[2]

############### get training data for the linear fits (skip if loading it)
num_ens_members = 90
projts = zeros((d, (L1+L2), num_ens_members)) 
# hist_mean_temp = zeros((1, (L1)*num_ens_members))
ens_gmt = zeros((num_ens_members, Int((L1+L2)/12)))

for i in ProgressBar(1:num_ens_members)
    files = ["data/historical/ts/ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc",
        # "data/RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
        # "data/RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
        "data/RCP45/ts/ts_Amon_MPI-ESM_rcp45_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
    tmps = ncData(files[1], "ts")
    projts[:, 1:L1, i] = project_timeseries(tmps.data, basis)
    ens_gmt[i, 1:Int(L1/12)] = get_gmt_list(tmps) 
    tmps = ncData(files[2], "ts")
    projts[:, L1+1:end, i] = project_timeseries(tmps.data, basis)
    ens_gmt[i, Int(L1/12)+1:end] = get_gmt_list(tmps)
end


# hfile = h5open("data/training_data_with45_20d_90ens.hdf5", "w")
# write(hfile, "projts", projts)
# write(hfile, "ens_gmt", ens_gmt)
# write(hfile, "num_ens_members", num_ens_members)
# close(hfile)