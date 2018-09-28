using Test
using Dorothy
using .Selectors
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Selection" begin

    m1 = readf("$(datapath)/1BTL.pdb")
    s1 = view(m1, Water)
    @test s1.ids == 2039:2237
    @test length(view(m1, Protein)) == 2032
    @test length(view(m1, Water)) == 199
    @test length(view(m1, Water | Protein)) == 2032 + 199
    @test length(view(m1, Water & Protein)) == 0
    @test length(view(m1, Expand(Index(1), by=Residue))) == 10
    @test length(view(m1, Expand(Index(1), by=Chain))) == 2236
    v1 = view(m1, Water & Within(5.0, Nter, pbc=false))
    @test length(v1) == 7
    v2 = view(m1, Water & Within(5.0, Nter, pbc=true))
    @test length(v2) == 7
    # m2 = readf("$(datapath)/PLC.pdb")
    # v2 = view(m2, Water & Within(3.0, Lipid, pbc=false))
    # writef("test.pdb", v2)

end # @testset
