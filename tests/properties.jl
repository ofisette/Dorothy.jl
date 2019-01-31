using Test
using LinearAlgebra
using Dorothy
using Dorothy.Graphs
using Dorothy.Geometry
using Dorothy.PBC
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Properties" begin

	@testset "Guessing elements" begin
		m1 = readf("$(datapath)/1BTL.gro")
		@test guesselement("CA", "CA") == "Ca"
		@test guesselement("CL", "CL") == "Cl"
		@test guesselement("CA", "ASP") == "C"
		@test guesselement("N", "ASP") == "N"
		@test_throws Exception guesselement("", "ASP")
		guesselements!(m1)
		@test m1.elements[1:4] == ["N", "C", "C", "O"]
	end

	@testset "Guessing masses" begin
		m1 = readf("$(datapath)/1BTL.gro")
		@test guessmass("MN1", "") == 0.0
		@test guessmass("HG11", "H") == 1.008
		@test_throws Exception guessmass("", "")
		guessmasses!(m1)
		@test m1.masses[1:4] == [14.00, 12.01, 12.01, 15.99]
	end

end # @testset
