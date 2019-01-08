using Test
using Dorothy
using Dorothy.Geometry
using Dorothy.PBC
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "..", "data")

@testset "PDB" begin

	@testset "Read" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		@test length(m1) == 2236
		sides, angles = pbcgeometry(m1.header.cell)
		@test sides ≈ [43.1, 64.4, 91.2]
		@test m1.header.title ==
				"CRYSTAL STRUCTURE OF ESCHERICHIA COLI TEM1 BETA-LACTAMASE" *
				" AT 1.8 ANGSTROMS RESOLUTION"
		@test m1.resids[1] == 26
		@test m1.resnames[1] == "HIS"
		@test m1.names[1] == "N"
		@test m1.ids[1] == 1
		@test m1.chainids[1] == "A"
		@test m1.occupancies[1] == 1.0
		@test m1.bfactors[1] == 14.51
		@test m1.R[end] ≈ [-12.176, -11.865, 48.715]
		@test ! haskey(m1, :V)
	end

	@testset "Write" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		io = IOBuffer()
		write(specify(io, "structure/x-pdb"), m1)
		seekstart(io)
		m2 = read(specify(io, "structure/x-pdb"))
		@test length(m2) == 2236
		sides, angles = pbcgeometry(m2.header.cell)
		@test sides ≈ [43.1, 64.4, 91.2]
		@test m2.header.title ==
				"CRYSTAL STRUCTURE OF ESCHERICHIA COLI TEM1 BETA-LACTAMASE" *
				" AT 1.8 ANGSTROMS RESOLUTION"
		@test m2.resids[1] == 26
		@test m2.resnames[1] == "HIS"
		@test m2.names[1] == "N"
		@test m2.ids[1] == 1
		@test m2.chainids[1] == "A"
		@test m2.occupancies[1] == 1.0
		@test m2.bfactors[1] == 14.51
		@test m2.R[end] ≈ [-12.176, -11.865, 48.715]
		@test ! haskey(m2, :V)
	end

	@testset "Convert" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		io = IOBuffer()
		write(specify(io, "structure/x-gro"), m1)
		seekstart(io)
		m2 = read(specify(io, "structure/x-gro"))
		@test length(m2) == 2236
		sides, angles = pbcgeometry(m2.header.cell)
		@test sides ≈ [43.1, 64.4, 91.2]
		@test m2.header.title ==
				"CRYSTAL STRUCTURE OF ESCHERICHIA COLI TEM1 BETA-LACTAMASE" *
				" AT 1.8 ANGSTROMS RESOLUTION"
		@test m2.resids[1] == 26
		@test m2.resnames[1] == "HIS"
		@test m2.names[1] == "N"
		@test m2.ids[1] == 1
		@test m2.R[end] ≈ [-12.18, -11.87, 48.72]
		@test ! haskey(m2, :V)
	end

end # @testset
