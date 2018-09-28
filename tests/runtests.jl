using Test
using Dorothy
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Dorothy" begin

    include("graphs.jl")
    include("geometry.jl")
    include("models.jl")
    include("guessing.jl")
    include("selection.jl")

    @testset "Formats" begin
        include("formats/gromos87.jl")
        include("formats/pdb.jl")
        include("formats/trr.jl")
        include("formats/xtc.jl")
    end

end # @testset
