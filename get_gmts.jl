using HDF5#, ProgressBars
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl")

L1, L2 = 1980, 1032 #for CMIP6

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]

for scenario in scenarios
    ############### get training data for the linear fits 
    num_ens_members = 50 # number of model runs used to train the emulator
    ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12))) #GMT in any case

    file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/"
    errors = []
    for i in 1:num_ens_members#ProgressBar(1:num_ens_members)
        try
            println("working on ensemble member $(i)")
            flush(stdout)
            files = file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"
            tmps = ncData(files, "tas")
            ens_gmt[i, :] = get_gmt_list(tmps) 
        catch
            println("missing values in ensemble member $(i)")
            push!(errors, i)
            flush(stdout)
        end
    end

    # for i in reverse(errors)
    # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run
    ens_gmt = ens_gmt[1:end .!= errors[:], :]
    # # end
    num_ens_members = size(ens_gmt)[1]

    hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "ens_gmt", ens_gmt)    
    write(hfile, "num_ens_members", num_ens_members)
    close(hfile)
end