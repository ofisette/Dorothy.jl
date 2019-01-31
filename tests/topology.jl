using Test
using LinearAlgebra
using Dorothy
using Dorothy.Graphs
using Dorothy.Geometry
using Dorothy.PBC
using Dorothy.Topology
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Topology" begin

	@testset "Automated assessment" begin
		m1 = readf("$(datapath)/1BTL.gro")
		topo = guesstopology!(m1)
		@test topo[1] == [2]
		@test topo[2] == [1,3,5]
		@test topo[3] == [2,4,11]
		@test topo[4] == [3]
		m2 = readf("$(datapath)/MHC.pdb")
		topo = guesstopology!(m2)
		@test topo[1] == []
		@test topo[3] == [4,5,6,7]
		@test topo[4] == [3]
		@test topo[6540] == [6541, 6542]
	end

end # @testset
