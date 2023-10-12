using Dates

function step_through(d::ncData, refs, assignment_func) #rght now it's just monthly
    #returns a dictionary, keyed on month, of markov chains for this dataset
    #refs should be of the format to fit the assignment function
    mcs = Dict()
    for t in 1:Int(length(d.timevec)/4) #CHANGE BACK TO FULL LENGTH
        month = Dates.month(d.timevec[t])
        class = assignment_func(d.data[:,:,t], refs[month])
        try
            push!(mcs[month], class)
        catch
            mcs[month] = [class]
        end
    end
    return mcs
end

##### assignment functions

#function to find closest image
function compare(img, refs)
    #returns index of closest ref 
    size(img) == size(refs[1]) || throw(ArgumentError("Dimension mismatch"))
    # let's use MSE
    errors = [mse(img, ref) for ref in refs]
    return argmin(errors)
end

#helper
function mse(img, ref)
    sq_diff = (img .- ref) .^ 2
    return mean(sq_diff)
end

######
#to do it via eof:
#our 'ref' here is a matrix U, columns are a basis for the data

function eof_assignment(snap, U) #this isn't what we want actually
    #snap is just a flat image
    v = reshape_data(snap)
    proj = projection(v, U)
    return argmax(proj)
end

#helper
function projection(v, U)
    dim = size(U,2)
    proj = zeros(dim)
    for i in 1:dim
        proj[i] = dot(v, U[:, i])/dot(U[:,i], U[:,i])
    end
    return proj
end
