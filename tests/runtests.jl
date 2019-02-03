using Test
using Dorothy
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Dorothy" begin

	include("graphs.jl")
	include("geometry.jl")
	include("pbc.jl")

	include("models.jl")
	include("properties.jl")
	include("hierarchies.jl")
	include("selections.jl")

	@testset "Formats" begin
		include("formats/gromos87.jl")
		include("formats/pdb.jl")
		include("formats/trr.jl")
		include("formats/xtc.jl")
		include("formats/ndx.jl")
	end

	include("topology.jl")
	include("ss.jl")

end # @testset
