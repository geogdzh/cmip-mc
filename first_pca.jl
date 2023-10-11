using NetCDF, Dates, Statistics, CairoMakie

include("utils.jl")

function flatten(d::ncData, start_year, end_year)
    #start through end year inclusive
    first = findfirst(x -> x == DateTime(start_year, 1, 1), d.timevec)
    last = findfirst(x -> x == DateTime(end_year, 12, 31), d.timevec)
    
end