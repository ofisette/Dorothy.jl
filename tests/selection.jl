using Test
using Dorothy
using .Selectors
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Selection" begin

	@testset "1BTL" begin # Small orthorhombic system
		model = readf("$(datapath)/1BTL.pdb")
		sel1 = view(model, Water)
		@test sel1.ids == 2039:2237
		@test length(view(model, Protein)) == 2032
		@test length(view(model, Water)) == 199
		@test length(view(model, Water | Protein)) == 2032 + 199
		@test length(view(model, Water & Protein)) == 0
		@test length(view(model, Expand(Index(1), by=Residue))) == 10
		@test length(view(model, Expand(Index(1), by=Chain))) == 2236
		sel2 = view(model, Water & Within(5.0, of=Nter))
		@test length(sel2) == 7
		# TODO: Increase the distance to get water from a periodic image, so we
		# can compare to the non-periodic system.
		sel3 = view(model, Water & Within(5.0, of=Nter))
		@test length(sel3) == 7
	end

	@testset "1BTL-np" begin # Non-periodic system
		model = readf("$(datapath)/1BTL.pdb")
		delete!(model.header, :cell)
	end

	@testset "MHC" begin # Rhombododecahedral system
		model = readf("$(datapath)/MHC.pdb")
	end

	@testset "PLC" begin # Very large orthorhombic system
		model = readf("$(datapath)/PLC.pdb")
		sel1 = view(model, Water & Expand(Within(3.0, of=Lipid), by=Residue))
	end

end # @testset
