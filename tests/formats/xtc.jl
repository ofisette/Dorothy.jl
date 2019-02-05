using Test
using Dorothy

@DorothyAll()

datapath = joinpath(@__DIR__, "..", "..", "data")

@testset "XTC" begin

	@testset "Single" begin
		f1 = infer("$(datapath)/eMHC.xtc")
		@test getformat(f1) == "trajectory/x-xtc"
		f2 = infer(open("$(datapath)/eMHC.xtc"))
		@test getformat(f2) == "trajectory/x-xtc"
		m1 = read(f2)
		@test m1.header.step == 0
		@test m1.header.time == 0.0
		@test length(m1) == 6002
		@test volume(m1.header.cell) ≈ 665416.28
		m2 = read!(f2, MolecularModel())
		@test m2.header.step == 500e3
		@test m2.header.time == 2e3
		@test isapprox(m2.R[1], [88.65, 82.36, 47.93], atol=0.01)
		@test isapprox(m2.R[end], [72.15, 56.0, 41.75], atol=0.01)
	end

	@testset "Trajectory" begin
		streamf("$(datapath)/eMHC.xtc") do s
			frame = MolecularModel()
			i = 0
			while !eof(s)
				read!(s, frame)
				i += 1
				if frame.header.step == 500e3
					@test frame.header.time == 2e3
					@test isapprox(frame.R[1], [88.65, 82.36, 47.93], atol=0.01)
					@test isapprox(frame.R[end], [72.15, 56.0, 41.75],
							atol=0.01)
				end
			end
			@test i == 501
		end
	end

end # @testset
