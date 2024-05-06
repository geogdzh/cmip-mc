using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

# L1, L2 = 1872, 1128 #for CMIP5
L1, L2 = 1980, 1032 #for CMIP6
scenario = "ssp585"
# scenario = "ssp119"

using_precip = true

######### load in basis
hfile = using_precip ? h5open("data/temp_precip/temp_precip_basis_1000d.hdf5", "r") : h5open("data/only_temp/temp_basis_1000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
close(hfile)
d = parse(Int, ARGS[1])
basis = basis[:, 1:d]

########## load in training data
hfile = using_precip ?  h5open("data/temp_precip/training_data_withpr_$(scenario)_1000d_49ens.hdf5", "r") : h5open("data/only_temp/training_data_$(scenario)_200d_49ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
ens_gmt = read(hfile, "ens_gmt")
num_ens_members = read(hfile, "num_ens_members")
close(hfile)

########## get the emulator itself
println("working on emulator for $(scenario) with $(d) dimensions")
flush(stdout)
ens_gmt = mean(ens_gmt, dims=1)
mean_coefs = get_mean_coefs(ens_projts, ens_gmt, degree=2)
# # covs = gmt_cov(ens_projts, ens_gmt)
# println("getting the covariances")
# flush(stdout)
# gmt_cov(ens_projts, ens_gmt, "$(scenario)_$(d)d")
# # chol_coefs = get_chol_coefs(covs, ens_gmt) #old version, which is faster
# println("getting the cholesky decomposition")
# flush(stdout)
# get_chols("$(scenario)_$(d)d")
println("done with the cholesky decomposition, now fitting the emulator")
flush(stdout)
chol_coefs = get_chol_coefs(ens_gmt, "$(scenario)_$(d)d")

hfile = using_precip ? h5open("data/temp_precip/gaussian_emulator_withpr_$(scenario)_$(d)d.hdf5", "w") : h5open("data/temp_only/gaussian_emulator_$(scenario)_$(d)d.hdf5", "w")
write(hfile, "mean_coefs", mean_coefs)
write(hfile, "chol_coefs", chol_coefs)
write(hfile, "basis", basis)
write(hfile, "num_ens_members", num_ens_members)
write(hfile, "ens_gmt", ens_gmt)
close(hfile)