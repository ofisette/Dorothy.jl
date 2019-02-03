using Test
using LinearAlgebra
using Dorothy
using Dorothy.Graphs
using Dorothy.Geometry
using Dorothy.PBC
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Hierarchies" begin

	@testset "Chains" begin
		m1 = MolecularModel(10)
		m1.chainids =
				[fill("A", 3)..., fill("B", 2)..., fill("C", 4)..., "D"]
		chainiter = chains(m1)
		@test [length(chain) for chain in chainiter] == [3,2,4,1]
		@test parentindices(chainiter[1]) == 1:3
		@test parentindices(chainiter[2]) == 4:5
		@test parentindices(chainiter[3]) == 6:9
		@test parentindices(chainiter[4]) == 10:10
		@test chainat(m1, 1) == view(m1, 1:3)
		@test chainat(m1, 2) == view(m1, 1:3)
		@test chainat(m1, 3) == view(m1, 1:3)
		@test chainat(m1, 4) == view(m1, 4:5)
		@test chainat(m1, 5) == view(m1, 4:5)
		@test chainat(m1, 6) == view(m1, 6:9)
		@test chainat(m1, 7) == view(m1, 6:9)
		@test chainat(m1, 8) == view(m1, 6:9)
		@test chainat(m1, 9) == view(m1, 6:9)
		@test chainat(m1, 10) == view(m1, 10:10)
	end

	@testset "Residues" begin
		m1 = MolecularModel(20)
		m1.resids = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5]
		m1.resnames = fill("WAT", 20)
		@test [length(res) for res in residues(m1)] == [4,4,4,4,4]
		m1.chainids = [fill("A", 10)..., fill("B", 10)...]
		resiter = residues(m1)
		@test [length(res) for res in resiter] == [4,4,2,2,4,4]
		@test parentindices(resiter[1]) == 1:4
		@test parentindices(resiter[2]) == 5:8
		@test parentindices(resiter[3]) == 9:10
		@test parentindices(resiter[4]) == 11:12
		@test parentindices(resiter[5]) == 13:16
		@test parentindices(resiter[6]) == 17:20
		@test residueat(m1, 1) == view(m1, 1:4)
		@test residueat(m1, 2) == view(m1, 1:4)
		@test residueat(m1, 3) == view(m1, 1:4)
		@test residueat(m1, 4) == view(m1, 1:4)
		@test residueat(m1, 5) == view(m1, 5:8)
		@test residueat(m1, 16) == view(m1, 13:16)
		@test residueat(m1, 17) == view(m1, 17:20)
		@test residueat(m1, 18) == view(m1, 17:20)
		@test residueat(m1, 19) == view(m1, 17:20)
		@test residueat(m1, 20) == view(m1, 17:20)
	end

	@testset "Fragments" begin
		m1 = MolecularModel(10)
		m1.topology = []
		pair!(m1.topology, (1, 2), (2, 3), (4, 5), (4, 6), (4, 7), (4, 8))
		fragiter = fragments(m1)
		@test [length(frag) for frag in fragiter] == [3,5,1,1]
		@test parentindices(fragiter[1]) == 1:3
		@test parentindices(fragiter[2]) == 4:8
		@test parentindices(fragiter[3]) == 9:9
		@test parentindices(fragiter[4]) == 10:10
		@test fragmentat(m1, 1) == view(m1, 1:3)
		@test fragmentat(m1, 2) == view(m1, 1:3)
		@test fragmentat(m1, 3) == view(m1, 1:3)
		@test fragmentat(m1, 4) == view(m1, 4:8)
		@test fragmentat(m1, 5) == view(m1, 4:8)
		@test fragmentat(m1, 6) == view(m1, 4:8)
		@test fragmentat(m1, 7) == view(m1, 4:8)
		@test fragmentat(m1, 8) == view(m1, 4:8)
		@test fragmentat(m1, 9) == view(m1, 9:9)
		@test fragmentat(m1, 10) == view(m1, 10:10)
	end

end # @testset
