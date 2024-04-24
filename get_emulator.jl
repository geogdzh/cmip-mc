using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

# L1, L2 = 1872, 1128 #for CMIP5
L1, L2 = 1980, 1032 #for CMIP6
scenario = "ssp585"
# scenario = "ssp119"

######### load in basis
hfile = h5open("data/temp_basis_1000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = 200
basis = basis[:, 1:d]

########## load in training data
hfile = h5open("data/training_data_$(scenario)_200d_49ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
ens_gmt = read(hfile, "ens_gmt")
num_ens_members = read(hfile, "num_ens_members")
close(hfile)

########## get the emulator itself
ens_gmt = mean(ens_gmt, dims=1)
mean_coefs = get_mean_coefs(ens_projts, ens_gmt)
covs = gmt_cov(ens_projts, ens_gmt)
chol_coefs = get_chol_coefs(covs, ens_gmt)

hfile = h5open("data/gaussian_emulator_$(scenario)_$(d)d.hdf5", "w")
write(hfile, "mean_coefs", mean_coefs)
write(hfile, "chol_coefs", chol_coefs)
write(hfile, "basis", basis)
write(hfile, "num_ens_members", num_ens_members)
write(hfile, "ens_gmt", ens_gmt)
close(hfile)