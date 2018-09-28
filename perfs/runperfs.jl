using Dorothy
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

include("trr.jl")
include("xtc.jl")
