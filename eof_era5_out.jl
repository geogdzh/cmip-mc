using NetCDF, Dates, Statistics, CairoMakie, LinearAlgebra, JLD
include("utils.jl")


twoyear = JLD.load("data/test_era5.jld", "goal_data_array")

