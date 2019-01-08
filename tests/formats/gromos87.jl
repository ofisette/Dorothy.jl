using Test
using Dorothy
using Dorothy.Geometry
using Dorothy.PBC
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "..", "data")

@testset "Gromos87" begin

	@testset "Read" begin
		m1 = readf("$(datapath)/1BTL.gro")
		@test length(m1) == 2236
		sides, angles = pbcgeometry(m1.header.cell)
		@test sides ≈ [43.1, 64.4, 91.2]
		@test m1.header.title == "BETA-LACTAMASE TEM1"
		@test m1.resids[1] == 26
		@test m1.resnames[1] == "HIS"
		@test m1.names[1] == "N"
		@test m1.ids[1] == 1
		@test m1.R[end] ≈ [-12.18, -11.86, 48.72]
		@test ! haskey(m1, :V)
	end

	@testset "Write" begin
		m1 = readf("$(datapath)/1BTL.gro")
		io = IOBuffer()
		write(specify(io, "structure/x-gro"), m1)
		seekstart(io)
		m2 = read(specify(io, "structure/x-gro"))
		@test length(m2) == 2236
		sides, angles = pbcgeometry(m2.header.cell)
		@test sides ≈ [43.1, 64.4, 91.2]
		@test m2.header.title == "BETA-LACTAMASE TEM1"
		@test m2.resids[1] == 26
		@test m2.resnames[1] == "HIS"
		@test m2.names[1] == "N"
		@test m2.ids[1] == 1
		@test m2.R[end] ≈ [-12.18, -11.86, 48.72]
		@test ! haskey(m2, :V)
	end

end # @testset
