using CairoMakie, ProgressBars, HDF5
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 


scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]
scenario = scenarios[1]
ens_ind = 8
for ens_ind in 1:50
    println(ens_ind)
    filehead = "/net/fs06/d3/mgeo/CMIP6/interim/"
    filetail = "$(scenario)/pr/r$(ens_ind)i1p1f1_$(scenario)_pr.nc"
    file = filehead * filetail

    pr = ncData(file, "pr")
    # tas.timevec
    # timeseries = month_to_year_avg(weighted_avg(tas))

    println(count_nan(pr.data)  )
end

lines(timeseries)
heatmap(tas.data[:,:,2])

function count_nan(arr::AbstractArray)
    count = 0
    locs = []
    for i in 1:M
        for j in 1:N
            for k in 1:L1
                value = arr[i,j,k]
                if ismissing(value)
                    count += 1
                    push!(locs, (i,j,k))
                end
            end
        end
    end
    return count, locs
end

function count_nan(arr::AbstractArray)
    count = 0
    for value in arr
        if ismissing(value)
            count += 1
        end
    end
    return count
end

cnt, locs = count_nan(tas.data) #221184


hfile = h5open("data/CMIP5/training_data_rcp85_20d_90ens.hdf5", "r")
projts_20 = read(hfile, "projts")
ens_gmt_20 = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/CMIP5/training_data_rcp85_100d_90ens.hdf5", "r")
projts_100 = read(hfile, "projts")
ens_gmt_100 = read(hfile, "ens_gmt")
close(hfile)

ens_gmt_20 == ens_gmt_100
projts_100[1:20,:,:] == projts_20

scenario = "ssp585"
ds = Dataset("/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/tas/r$(1)i1p1f1_$(scenario)_tas.nc")
tas = ds["tas"][:]
close(ds)

tas = ncData("/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/tas/r$(10)i1p1f1_$(scenario)_tas.nc", "tas")
gmt_119 = get_gmt_list(tas)
gmt_585 = get_gmt_list(tas)
lines(gmt_119)
lines(gmt_585)

#

file = "/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/r1i1p1f1/historical/pr/250_km/mon/1850/CMIP6_MPI-ESM1-2-LR_r1i1p1f1_historical_pr_250_km_mon_gn_1850.nc"
pr = ncData(file, "pr")

using NCDatasets
ds = Dataset(file)

ncgen(file, "samplepr.jl")

hfile = h5open("data/temp_basis_1000d.hdf5")
temp_basis = read(hfile, "basis")
close(hfile)

hfile = h5open("data/temp_precip_basis_1000d.hdf5")
temp_precip_basis = read(hfile, "basis")
close(hfile)

temp_basis == temp_precip_basis[1:M*N, :]


#ok, so it has to be trained together


hfile = h5open("data/temp_precip/ens_vars_withpr_ssp119.hdf5", "r")
ens_means_pr_20 = read(hfile, "ens_means_pr_20")
ens_means_pr_100 = read(hfile, "ens_means_pr_100")  
ens_vars_pr_20 = read(hfile, "ens_vars_pr_20")
ens_vars_pr_100 = read(hfile, "ens_vars_pr_100")
close(hfile)

hfile = h5open("data/temp_precip/ens_vars_withpr_ssp585.hdf5", "r")
b_ens_means_pr_20 = read(hfile, "ens_means_pr_20")
b_ens_means_pr_100 = read(hfile, "ens_means_pr_100")  
b_ens_vars_pr_20 = read(hfile, "ens_vars_pr_20")
b_ens_vars_pr_100 = read(hfile, "ens_vars_pr_100")
close(hfile)

heatmap(b_ens_vars_pr_20[:,:,401])
ens_vars_pr_20[:,:,4] == b_ens_vars_pr_20[:,:,4]


#
hfile = h5open("data/temp_precip/gaussian_emulator_withpr_ssp585_200d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
close(hfile)


# hfile = h5open("data/temp_only/vars_$(scenario)_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
# true_var_temp = read(hfile, "true_var")
# num_ens_members = read(hfile, "num_ens_members")
# true_ens_gmt = read(hfile, "true_ens_gmt")
# true_ens_mean_temp = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
# close(hfile)

# hfile = h5open("data/temp_precip/vars_pr_$(scenario)_50ens.hdf5", "r") #this is the actual ensemble variance of the CMIP model
# true_var_pr = read(hfile, "true_var")
# # num_ens_members = read(hfile, "num_ens_members")
# # true_ens_gmt = read(hfile, "true_ens_gmt")
# true_ens_mean_pr = read(hfile, "true_ens_mean")[:,:,:,1] #NEW!
# close(hfile)


# hfile = h5open("data/temp_precip/ens_vars_withpr_$(scenario).hdf5", "r")
# ens_means_pr_20 = read(hfile, "ens_means_pr_20")
# ens_means_pr_100 = read(hfile, "ens_means_pr_100")  
# ens_vars_pr_20 = read(hfile, "ens_vars_pr_20")
# ens_vars_pr_100 = read(hfile, "ens_vars_pr_100")
# ens_means_temp_20 = read(hfile, "ens_means_tas_20")
# ens_means_temp_100 = read(hfile, "ens_means_tas_100")  
# ens_vars_temp_20 = read(hfile, "ens_vars_tas_20")
# ens_vars_temp_100 = read(hfile, "ens_vars_tas_100")
# close(hfile)


hfile = h5open("data/temp_precip/training_data_withpr_ssp585_20d_49ens.hdf5")
projts_20 = read(hfile, "projts")
ens_gmt_20 = read(hfile, "ens_gmt")
close(hfile)

size(projts_20)
projts_20[:,:,1]

hfile = h5open("data/process/test_emulator_withpr_ssp585_20d.hdf5", "r") #this is to test the new chol/cov offloading
# mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
close(hfile)

hfile = h5open("data/temp_precip/gaussian_emulator_withpr_ssp585_20d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs")
old_chol_coefs = read(hfile, "chol_coefs")
old_basis = read(hfile, "basis")
close(hfile)

chol_coefs == old_chol_coefs
basis == old_basis

# looks like new method works just fine!

hfile = h5open("data/process/chols_ssp585_1000d.hdf5", "r")
chols = read(hfile, "chols_251")
D = read(hfile, "D")
num_years  = read(hfile, "num_years")
close(hfile)