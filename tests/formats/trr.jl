using Test
using Dorothy
using ..Dorothy.Geometry
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "..", "data")

@testset "TRR" begin

    @testset "Single" begin
        f1 = guess("$(datapath)/eMHC.trr")
        @test getformat(f1) == "trajectory/x-trr"
        f2 = guess(open("$(datapath)/eMHC.trr"))
        @test getformat(f2) == "trajectory/x-trr"
        m1 = read(f2)
        @test m1.header.step == 0
        @test m1.header.time == 0.0
        @test length(m1) == 6002
        @test volume(m1.header.cell) ≈ 665416.28
        m2 = read!(f2, MolecularModel())
        @test m2.header.step == 500e3
        @test m2.header.time == 2e3
        @test isapprox(m2.R[:,1], [88.65, 82.36, 47.93], atol=0.01)
        @test isapprox(m2.R[:,end], [72.15, 56.0, 41.75], atol=0.01)
    end

    @testset "Trajectory" begin
        streamf("$(datapath)/eMHC.trr") do s
            @test length(s) == 501
            frame = MolecularModel()
            i = 0
            while !eof(s)
                read!(s, frame)
                i += 1
                if frame.header.step == 500e3
                    @test frame.header.time == 2e3
                    @test isapprox(frame.R[:,1], [88.65, 82.36, 47.93],
                            atol=0.01)
                    @test isapprox(frame.R[:,end], [72.15, 56.0, 41.75],
                            atol=0.01)
                end
            end
            @test i == length(s)
        end
    end

    @testset "Round-trip" begin
        buffer1 = IOBuffer(read("$(datapath)/eMHC.trr"))
        buffer2 = IOBuffer(write=true)
        traj1 = streamf(buffer1)
        traj2 = streamf(specify(buffer2, "trajectory/x-trr"),
                nparticles=traj1.meta.nparticles)
        frame = MolecularModel()
        while !eof(traj1)
            read!(traj1, frame)
            write(traj2, frame)
        end
        @test take!(buffer1) == take!(buffer2)
    end

end # @testset
