using Test
using Dorothy
using Dorothy.Hierarchies
using Dorothy.PBC
using .Selectors
using Formats
using FormatStreams
using FormatCodecs

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
		@test length(view(model, Index(1, by=Residue))) == 10
		@test length(view(model, Expand(Index(1), by=Residue))) == 10
		@test length(view(model, Index(1, by=Chain))) == 2236
		@test length(view(model, Expand(Index(1), by=Chain))) == 2236
		sel2 = view(model, Water & Within(5.0, of=Nter))
		@test length(sel2) == 7
		sel3 = view(model, Water & Within(45.0, of=Nter))
		@test length(sel3) == 199
	end

	@testset "1BTL-np" begin # Non-periodic system
		model = readf("$(datapath)/1BTL.pdb")
		delete!(model.header, :cell)
		sel1 = view(model, Water & Within(45.0, of=Nter))
		@test length(sel1) == 188
	end

	@testset "MHC" begin # Rhombododecahedral system
		model = readf("$(datapath)/MHC.pdb")
		sel1 = view(model, Water & Within(5.0, of=Protein))
		@test length(sel1) == 6475
		sel2 = view(model, Cter)
		@test length(sel2) == 56
		@test length(residues(sel2)) == 3
		sel3 = view(model, Nter)
		@test length(sel3) == 57
		@test length(residues(sel3)) == 3
		bitvec = falses(length(model))
		for i = 1:6539
			bitvec[i] = true
		end
		sel4 = view(model, Map(bitvec))
		@test length(sel4) == 6539
		sel5 = view(model, Protein)
		@test sel5 == sel4
		sel6 = view(model, Index(1:11))
		@test length(sel6) == 11
		@test all(sel6.resnames .== "GLY")
		sel7 = view(model, Index(1))
		@test sel7[].name == "MN1"
		sel8 = view(model, Index(277, by=Residue))
		@test length(sel8) == 27
		@test all(sel8.resnames .== "ILE")
		sel9 = view(model, Index(3, by=Chain))
		@test length(sel9) == 149
		@test all(sel9.chainids .== "C")
		sel10 = view(model, Index(first, ineach=Chain))
		@test length(sel10) == 6
		Dorothy.Topology.infertopology!(model)
		@test length(fragments(view(model, 1:6539))) == 3
		sel11 = view(model, Index(3, by=Fragment))
		@test sel11 == sel9
		@test length(fragments(view(model, Protein))) == 3
		# Test expand and restrict, comparing with existing sels.
	end

	@testset "PLC" begin # Very large orthorhombic system
		model = readf("$(datapath)/PLC.pdb.bz2")
		sel1 = view(model, Water & Expand(Within(3.0, of=Lipid), by=Residue))
		@test length(sel1) == 30474
		sel2 = view(model, Water & Restrict(Within(3.0, of=Lipid), by=Residue))
		@test length(sel2) == 14271
	end

end # @testset
