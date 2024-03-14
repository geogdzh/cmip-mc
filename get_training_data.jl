using HDF5, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

L1 = 1872 #for historical (generalize!)
L2 = 1128 #for rcp (generalize!)

############## load a basis
# hfile = h5open("data/temp-precip_basis.hdf5", "r")
hfile = h5open("data/temp_basis_1000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = 600 #size(basis)[2]
basis = basis[:, 1:d]

############### get training data for the linear fits 
num_ens_members = 90 # number of model runs used to train the emulator
projts = zeros((d, (L1+L2), num_ens_members)) 
# hist_mean_temp = zeros((1, (L1)*num_ens_members))
ens_gmt = zeros((num_ens_members, Int((L1+L2)/12)))

file_head = "/net/fs06/d3/CMIP5/MPI-GE/"
for i in ProgressBar(1:num_ens_members)
    files = [file_head*"historical/ts/ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc",
        # file_head*"RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
        file_head*"RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
        # file_head*"RCP45/ts/ts_Amon_MPI-ESM_rcp45_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
    tmps = ncData(files[1], "ts")
    projts[:, 1:L1, i] = project_timeseries(tmps.data, basis)
    ens_gmt[i, 1:Int(L1/12)] = get_gmt_list(tmps) 
    tmps = ncData(files[2], "ts")
    projts[:, L1+1:end, i] = project_timeseries(tmps.data, basis)
    ens_gmt[i, Int(L1/12)+1:end] = get_gmt_list(tmps)
end


hfile = h5open("data/training_data_rcp85_$(d)d_90ens.hdf5", "w")
write(hfile, "projts", projts)
write(hfile, "ens_gmt", ens_gmt)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)

########## get the emulator itself
ens_gmt = mean(ens_gmt, dims=1)
ens_projts = projts
mean_coefs = get_mean_coefs(ens_projts, ens_gmt)
covs = gmt_cov(ens_projts, ens_gmt)
chol_coefs = get_chol_coefs(covs, ens_gmt)

hfile = h5open("data/gaussian_emulator_rcp85_$(d)d.hdf5", "w")
write(hfile, "mean_coefs", mean_coefs)
write(hfile, "chol_coefs", chol_coefs)
write(hfile, "basis", basis)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)