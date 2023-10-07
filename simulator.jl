using Dates

function step_through(d::ncData, refs) #rght now it's just monthly
    #returns a dictionary, keyed on month, of markov chains for this dataset
    mcs = Dict()
    for t in 1:Int(length(d.timevec)/4) #CHANGE BACK TO FULL LENGTH
        month = Dates.month(d.timevec[t])
        class = compare(d.data[:,:,t], refs[month])
        try
            push!(mcs[month], class)
        catch
            mcs[month] = [class]
        end
    end
    return mcs
end


#function to find closest image
function compare(img, refs)
    #returns index of closest ref 
    size(img) == size(refs[1]) || throw(ArgumentError("Dimension mismatch"))
    # let's use MSE
    errors = [mse(img, ref) for ref in refs]
    return argmin(errors)
end

function mse(img, ref)
    sq_diff = (img .- ref) .^ 2
    return mean(sq_diff)
end
