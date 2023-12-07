### training functions
# gmt = 290.
# vars1 = var_coefs[1, :, 2].*gmt .+ var_coefs[1, :, 1]
# vars2 = var_coefs[1+1, :, 2].*gmt .+ var_coefs[1+1, :, 1]
# vars = vcat(vars1, vars2)
# diagm(sqrt.(vars))*corrs[:,:,1]*diagm(sqrt.(vars))

function get_cov(gmt, corrs, var_coefs) #need all the corrs together
    # size(corrs) = (2d, 2d, 12)
    covs = zeros(size(corrs))
    for i in 1:12
        vars1 = var_coefs[i, :, 2].*gmt .+ var_coefs[i, :, 1] #first month
        vars2 = []
        try
            vars2 = var_coefs[i+1, :, 2].*gmt .+ var_coefs[i+1, :, 1] #second month
        catch
            vars2 = var_coefs[1, :, 2].*gmt .+ var_coefs[1, :, 1] #second month
        end
        vars = vcat(vars1, vars2)
        covs[:,:,i] = diagm(sqrt.(vars))*corrs[:,:,i]*diagm(sqrt.(vars))
    end
    return covs
end

function get_means(gmt, mean_coefs)
    out = zeros((size(mean_coefs)[2],12))
    for i in 1:12
        out[:,i] = mean_coefs[i,:,2].*gmt .+ mean_coefs[i, :, 1]
    end
    return out
end

### testing functions

function emulate(gmt_list, mean_coefs, corrs, var_coefs)
    num_years = length(gmt_list)
    # gmt assumed constant for each year
    trajectory = zeros(size(mean_coefs)[2], 12*num_years) 

    gmt = gmt_list[1]
    covs = get_cov(gmt, corrs, var_coefs)
    means = get_means(gmt, mean_coefs) #list of twelve

    dec = means[:,12]
    trajectory[:,1] = emulate_step(dec, 12, covs, means; new_means=means) #think again aobut whether or not this makes sense

    for year in 1:num_years #year is an index
        for i in 1:11 #here i is the prev_month
            trajectory[:, (i+1)+(year-1)*12] = emulate_step(trajectory[:, i+(year-1)*12], i, covs, means)
        end
        if year != num_years #if that wasn't the last year
            #now set up for the next year
            gmt = gmt_list[year+1]
            new_means = get_means(gmt, mean_coefs)
            trajectory[:, (12+1)+(year-1)*12] = emulate_step(trajectory[:, 12+(year-1)*12], 12, covs, means; new_means=new_means) #do the december transition
            covs = get_cov(gmt, corrs, var_coefs)
            means = new_means
        end
    end
    return trajectory
end


function emulate_step(prev_val, prev_month, covs, means; new_means=nothing)
    # prev_val - value of preceding month
    # prev_month - integer of preceding month
    # covs - matrix of covariances for the GMT of the preceding month; size = (10, 10, 12) where 10=2d
    # means - vector of means for the GMT of the preceding month
    # new_means - vector of means for the GMT of the new month, if different
    μ_1 = means[:,prev_month]
    μ_2 = isnothing(new_means) ? means[:,prev_month+1] : new_means[:,1] #hardcoded for switch to happen in january!
    Σ = covs[:,:,prev_month]
    d = Int(size(covs)[1]/2)
    Σ_11, Σ_12, Σ_21, Σ_22 = Σ[1:d,1:d], Σ[1:d,d+1:end], Σ[d+1:end,1:d], Σ[d+1:end,d+1:end]
    μ_out = μ_2 .+ Σ_21 * inv(Σ_11) * (prev_val .- μ_1)
    Σ_out = Σ_22 .- Σ_21 * inv(Σ_11) * Σ_12
    Σ_out = round.(Σ_out,digits=4)
    dist = MvNormal(μ_out, Σ_out)
    return rand(dist)
end