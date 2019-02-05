using Test
using Dorothy

@DorothyAll()

datapath = joinpath(@__DIR__, "..", "..", "data")

@testset "NDX" begin

	@testset "Read" begin
		index = readf("$(datapath)/PLC.ndx")
		@test index["system"] == 1:1612871
		@test index["protein"] == 1:81923
		@test index["Ahc_mod1"] == 1:5854
		@test index["B2m_mod1"] == 5855:7587
	end

	@testset "Write" begin
		index1 = readf("$(datapath)/PLC.ndx")
		io = IOBuffer()
		write(specify(io, "text/x-ndx"), index1)
		seekstart(io)
		index2 = read(specify(io, "text/x-ndx"))
		@test index2 == index1
	end

end # @testset
