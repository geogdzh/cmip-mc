## load the emulator
hfile = h5open("data/temp_precip/gaussian_emulator_withpr_ssp585_200d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
num_ens_members = read(hfile, "num_ens_members")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/temp_precip/training_data_withpr_ssp585_200d_49ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
# ens_gmt = read(hfile, "ens_gmt")
# num_ens_members = read(hfile, "num_ens_members")
close(hfile)


# hfile = h5open("data/temp_precip/gaussian_emulator_withpr_ssp119_20d.hdf5", "r")
# mean_coefs_119 = read(hfile, "mean_coefs")
# # chol_coefs_119 = read(hfile, "chol_coefs")
# # basis_119 = read(hfile, "basis")
# num_ens_members_119 = read(hfile, "num_ens_members")
# ens_gmt_119 = read(hfile, "ens_gmt")
# close(hfile)

# hfile = h5open("data/temp_precip/training_data_withpr_ssp119_20d_49ens.hdf5", "r")
# ens_projts_119 = read(hfile, "projts")[1:d, :, :]
# # ens_gmt = read(hfile, "ens_gmt")
# # num_ens_members = read(hfile, "num_ens_members")
# close(hfile)

hfile = h5open("data/ssp119_gmts_50ens.hdf5", "r")
ens_gmt_119 = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/temp_precip/projts_withpr_ssp119_200d_50ens.hdf5", "r")
projts_119 = read(hfile, "projts")
close(hfile)

begin
    fig = Figure(resolution=(2000, 2000))
    ax = Axis(fig[1,1])
    lines!(ax, ens_gmt[:])
    # lines!(ax, ens_gmt_119[:])
    display(fig)
end

####test
num_buckets = 10
bit = Int(size(ens_projts)[2]/num_buckets)
bitcorrs = zeros((2*d, 2*d, 12, num_buckets))
for i in 1:num_buckets
    bitcorrs[:,:,:,i] = get_corrs(ens_projts[:,(i-1)*bit+1:i*bit,:])
end
begin
    fig = Figure(resolution=(2000, 2000))
    month = 1
    ax = Axis(fig[1,1])
    heatmap!(ax, bitcorrs[d+1:end,1:d,month,5].-bitcorrs[d+1:end,1:d,month,1], colorrange=(-1,1),colormap=:balance)
    display(fig)
end
begin
    fig = Figure(resolution=(2000, 2000))
    # month = 1
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i,j], title = "mode $i vs mode $j of successive month")
            for month in 1:12
                lines!(ax, [x for x in 1:num_buckets], [bitcorrs[d+i,j,month,x] for x in 1:num_buckets])
            end
        end
    end
    # save("figs/corrs_over_time_10_buckets.png", fig)
    display(fig)
end

#mean fits
begin
    fig = Figure(resolution=(2000, 1000))
    # ax = Axis(fig[1,1])
    mon = 1
    # mode = 100
    for mode in 1:10
        # ax = Axis(fig[mode,1])
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[:], ens_projts[mode, mon:12:end, i], markersize=7, alpha=0.5, color=:red) 
            scatter!(ax, ens_gmt_119[:], ens_projts_119[mode, mon:12:end, i], markersize=7, alpha=0.5, color=:blue) 
        end
        ens_mean = mean(ens_projts[mode, mon:12:end, :], dims=2)
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], [mean_coefs[mon, mode, 2].*x.+mean_coefs[mon, mode, 1] for x in ens_gmt[:]], color=:black, linewidth=4)

        # ens_mean_119 = mean(ens_projts_119[mode, mon:12:end, :], dims=2)
        # scatter!(ax, ens_gmt_119[:], ens_mean_119[:], color=:brown, markersize=7)
        # lines!(ax, ens_gmt_119[:], [mean_coefs_119[mon, mode, 2].*x.+mean_coefs_119[mon, mode, 1] for x in ens_gmt_119[:]], color=:brown, linewidth=4)
    end
    # save("figs/jan_mean_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 

#var fits
#needs: vars, num_ens_members, var_coefs :: after an emulator has been trained

var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=true)
var_coefs_119, vars_119 = get_var_coefs(ens_projts_119, ens_gmt_119, mean_coefs_119; return_vars=true)


begin
    fig = Figure(resolution=(2000, 1000))
    # ax = Axis(fig[1,1])
    mon = 1
    # mode = 100
    for mode in 1:10
        # ax = Axis(fig[mode,1])
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="std for mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[:], sqrt.(vars[mon,mode,i,:]), markersize=7, alpha=0.5, color=:red)
            scatter!(ax, ens_gmt_119[:], sqrt.(vars_119[mon,mode,i,:]), markersize=7, alpha=0.5, color=:blue)
            # scatter!(ax, ens_gmt[:], (vars[mon,mode,i,:]), markersize=7, alpha=0.5) 
        end
        ens_mean = sqrt.(mean(dropdims(vars[mon:12:end, mode, :, :], dims=1), dims=1))
        # ens_mean = (mean(dropdims(vars[mon:12:end, mode, :, :], dims=1), dims=1))
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], sqrt.([var_coefs[mon, mode, 2].*x.+var_coefs[mon, mode, 1] for x in ens_gmt[:]]), color=:black, linewidth=4)
        # lines!(ax, ens_gmt[:], ([var_coefs[mon, mode, 2].*x.+var_coefs[mon, mode, 1] for x in ens_gmt[:]]), color=:black, linewidth=4)
    end
    # save("figs/jan_var_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 

#cov fits
begin
    fig = Figure(resolution=(1400, 1400))
    mon = 1
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i,j])
            x = hcat(fill(1., length(ens_gmt)), ens_gmt[:])
            y = covs[d+i,j,mon,:] #CHANGE HERE for corrs/covs # also it's lookin at the non-symmetrical part: covariance btw modes of diff months
            scatter!(ax, ens_gmt[:], y, color=:orange, alpha=0.5) 
            fits = x \ y
            # print(typeof(fit))
            lines!(ax, ens_gmt[:], fits[2].*ens_gmt[:] .+ fits[1], color=:black, linewidth=3)
            # lines!(ax, ens_gmt[:], [chol_coefs[month, i, j, 2].*x.+chol_coefs[month, i, j, 1] for x in ens_gmt[:]], color=:black, linewidth=3)
        end
    end
    save("figs/jan_cov_fits_ssp585.png", fig)
    display(fig)
end

#chol fits
begin
    fig = Figure(resolution=(1000, 1000))
    month = 1
    for i in 1:10
        for j in 1:10
            ax = Axis(fig[i,j])
            scatter!(ax, ens_gmt[:], chols[i,j,month,:], color=:orange, alpha=0.5) 
            lines!(ax, ens_gmt[:], [chol_coefs[month, i, j, 2].*x.+chol_coefs[month, i, j, 1] for x in ens_gmt[:]], color=:black, linewidth=3)
        end
    end
    display(fig)
end

