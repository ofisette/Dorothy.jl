using Test
using LinearAlgebra
using Dorothy
using Dorothy.Graphs
using Dorothy.Geometry
using Dorothy.PBC
using Dorothy.Properties
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Properties" begin

	@testset "Guessing elements" begin
		m1 = readf("$(datapath)/1BTL.gro")
		@test inferelement("CA", "CA") == "Ca"
		@test inferelement("CL", "CL") == "Cl"
		@test inferelement("CA", "ASP") == "C"
		@test inferelement("N", "ASP") == "N"
		@test_throws Exception inferelement("", "ASP")
		inferelements!(m1)
		@test m1.elements[1:4] == ["N", "C", "C", "O"]
	end

	@testset "Guessing masses" begin
		m1 = readf("$(datapath)/1BTL.gro")
		@test infermass("MN1", "") == 0.0
		@test infermass("HG11", "H") == 1.008
		@test_throws Exception infermass("", "")
		infermasses!(m1)
		@test m1.masses[1:4] == [14.00, 12.01, 12.01, 15.99]
	end

end # @testset
