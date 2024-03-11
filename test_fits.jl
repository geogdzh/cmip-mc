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
            scatter!(ax, ens_gmt[:], ens_projts[mode, mon:12:end, i], markersize=7, alpha=0.5) 
        end
        ens_mean = mean(ens_projts[mode, mon:12:end, :], dims=2)
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], [mean_coefs[mon, mode, 2].*x.+mean_coefs[mon, mode, 1] for x in ens_gmt[:]], color=:black, linewidth=4)
    end
    save("figs/jan_mean_fits_10_modes_rcp45.png", fig)
    display(fig)
end 

#var fits
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
            scatter!(ax, ens_gmt[:], sqrt.(vars[mon,mode,i,:]), markersize=7, alpha=0.5) 
        end
        ens_mean = sqrt.(mean(dropdims(vars[mon:12:end, mode, :, :], dims=1), dims=1))
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], sqrt.([var_coefs[mon, mode, 2].*x.+var_coefs[mon, mode, 1] for x in ens_gmt[:]]), color=:black, linewidth=4)
    end
    save("figs/jan_var_fits_10_modes_rcp45.png", fig)
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
            y = covs[d+i,j,mon,:] #CHANGE HERE for corrs/covs
            scatter!(ax, ens_gmt[:], y, color=:orange, alpha=0.5) 
            fits = x \ y
            # print(typeof(fit))
            lines!(ax, ens_gmt[:], fits[2].*ens_gmt[:] .+ fits[1], color=:black, linewidth=3)
            # lines!(ax, ens_gmt[:], [chol_coefs[month, i, j, 2].*x.+chol_coefs[month, i, j, 1] for x in ens_gmt[:]], color=:black, linewidth=3)
        end
    end
    save("figs/jan_cov_fits_rcp45.png", fig)
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

