function get_tvec(filename)
    #returns vector marking the first day of each month 
    # return DateTime(1850,1,1)+Day.(ncread(filename,"time"))
   [ DateTime(1850,1,1)+Month(x) for x in range(0,length(ncread(filename, "time"))-1)]
end

# function index(vec, value) 

#     return findfirst(x -> x==value, vec)
# end

function find_closest(vec, value)
    #finds index of closest value in vector
    geq = searchsortedfirst(vec, value)
    leq = searchsortedlast(vec, value)
    if abs(vec[leq]-value) < abs(vec[geq]-value)
        return leq
    else 
        return geq
    end
end

function convert_lon(lon)
    # in: lon on a +/- scale
    # out: lon on 360
    if lon > 0
        return lon
    else
        return 360 + lon
    end
end

function slice_map(d::ncData, lon1, lon2, lat1, lat2)
    lon1 = find_closest(d.lonvec, convert_lon(lon1))
    lon2 = find_closest(d.lonvec, convert_lon(lon2))
    lat1 = find_closest(d.latvec, lat1)
    lat2 = find_closest(d.latvec, lat2)
    if lon1 > lon2
        #we're passing over the prime meridian
        dataslice = vcat(d.data[lon1:end, lat1:lat2, :], d.data[1:lon2, lat1:lat2, :])
        newlonvec = vcat(d.lonvec[lon1:end], d.lonvec[1:lon2])
    else
        dataslice = d.data[lon1:lon2, lat1:lat2, :]
        newlonvec = d.lonvec[lon1:lon2]
    end
    newlatvec = d.latvec[lat1:lat2]
    return ncData(dataslice, newlonvec, newlatvec, d.timevec)
end

# Base.getindex(m::ncData, ind::Int) = slice #see abt implementing this later

struct ncData{D, NV, TV, T}
    data::D
    lonvec::NV
    latvec::TV
    timevec::T
end

function ncData(file, varname)
    data = ncread(file, varname)
    lonvec = ncread(file, "lon")
    latvec = ncread(file, "lat")
    timevec = [DateTime(1850,1,1)+Month(x) for x in range(0,length(ncread(file, "time"))-1)]
    return ncData(data, lonvec, latvec, timevec)
end

