using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

# L1, L2 = 1872, 1128 #for CMIP5
L1, L2 = 1980, 1032 #for CMIP6

############## load a basis
# hfile = h5open("data/temp-precip_basis.hdf5", "r")
hfile = h5open("data/temp_basis_1000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = 200 #size(basis)[2]
basis = basis[:, 1:d]

# scenario = "ssp119"
scenario = "ssp585"

############### get training data for the linear fits 
num_ens_members = 50 # number of model runs used to train the emulator
projts = zeros((d, (L1+L2), num_ens_members)) 
# hist_mean_temp = zeros((1, (L1)*num_ens_members))
ens_gmt = zeros((num_ens_members, Int((L1+L2)/12)))

# file_head = "/net/fs06/d3/CMIP5/MPI-GE/"
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
errors = []
for i in 1:num_ens_members#ProgressBar(1:num_ens_members)
    try
        println("working on ensemble member $(i)")
        flush(stdout)
        # files = [file_head*"historical/ts/ts_Amon_MPI-ESM_historical_r$(string(i, pad=3))i1850p3_185001-200512.nc",
            # file_head*"RCP26/ts/ts_Amon_MPI-ESM_rcp26_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
            # file_head*"RCP85/ts/ts_Amon_MPI-ESM_rcp85_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
            # file_head*"RCP45/ts/ts_Amon_MPI-ESM_rcp45_r$(string(i, pad=3))i2005p3_200601-209912.nc"]
        files = [ file_head*"historical/tas/r$(i)i1p1f1_historical_tas.nc",
            file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"]
        tmps = ncData(files[1], "tas")
        projts[:, 1:L1, i] = project_timeseries(tmps.data, basis)
        ens_gmt[i, 1:Int(L1/12)] = get_gmt_list(tmps) 
        tmps = ncData(files[2], "tas")
        projts[:, L1+1:end, i] = project_timeseries(tmps.data, basis)
        ens_gmt[i, Int(L1/12)+1:end] = get_gmt_list(tmps)
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

hfile = h5open("data/training_data_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
write(hfile, "projts", projts)
write(hfile, "ens_gmt", ens_gmt)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)
