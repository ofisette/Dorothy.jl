using Test
using LinearAlgebra
using Dorothy
using Dorothy.Graphs
using Dorothy.Geometry
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Models" begin

    @testset "Header" begin

        h1 = MolecularModelHeader()
        @test isempty(keys(h1))
        @test_throws Exception h1.smurf = true
        h1.time = 10.0
        @test h1.time == 10.0
        @test get(h1, :time, nothing) == 10.0
        @test length(keys(h1)) == 1
        h1.time = 10
        @test h1.time == 10.0
        @test length(keys(h1)) == 1
        @test_throws Exception h1.time = "ten dot zero"
        @test length(keys(h1)) == 1
        delete!(h1, :time)
        @test isempty(keys(h1))
        h1.title = "something"
        @test_throws Exception h1.title = nothing
        @test length(keys(h1)) == 1
        h1.step = 5
        @test_throws Exception h1.step = 2.5
        @test length(keys(h1)) == 2
        h1.lambda = 0.6
        @test length(keys(h1)) == 3
        h1.cell = pbccell(I)
        h1.cell = pbccell(25.0)
        @test sum(h1.cell, dims=1) ≈ [25.0 25.0 25.0]
        @test_throws Exception h1.cell = zeros(3,4)
        h1.virial = pbccell(I)
        h1.pressure = pbccell(I)
        show(IOBuffer(), h1)
        show(IOBuffer(), MIME"text/plain"(), h1)
        @test h1.lambda == h1[:lambda]
        @test !haskey(h1, :smurf)
        @test get(h1, :lambda, nothing) == 0.6
        @test get(h1, :smurf, nothing) == nothing
        @test get(() -> nothing, h1, :smurf) == nothing
        delete!(h1, :time)
        @test get!(h1, :time, 5.0) == 5.0
        @test h1.time == 5.0
        delete!(h1, :lambda)
        @test get!(() -> 1.0, h1, :lambda) == 1.0
        @test h1.lambda == 1.0

    end

    @testset "MolecularModel" begin

        m1 = MolecularModel()
        m2 = MolecularModel(10)
        @test isempty(m1)
        @test length(m2) == 10
        @test isempty(keys(m1))
        @test isempty(keys(m2))
        show(IOBuffer(), m2)
        show(IOBuffer(), MIME"text/plain"(), m2)
        @test firstindex(m2) == 1
        @test lastindex(m2) == 10
        @test collect(eachindex(m2)) == collect(1:10)
        @test_throws Exception m2.charges = fill("smurf", 10)
        m2.charges = fill(10.0, 10)
        @test m2.charges[1] == 10.0
        @test m2[1].charge == 10.0
        @test length(keys(m2)) == 1
        @test_throws Exception m2.occupancies = fill("tapir", 10)
        m2.occupancies = fill(true, 10)
        @test_throws Exception m2.bfactors = fill(nothing, 10)
        m2.bfactors = fill(10, 10)
        @test length(keys(m2)) == 3
        delete!(m2, :bfactors)
        @test length(keys(m2)) == 2
        delete!(m2, :charges)
        delete!(m2, :occupancies)
        @test isempty(keys(m2))
        @test_throws Exception m2.R = fill(0.0, (10,3))
        @test_throws Exception m2.R = fill(nothing, (3,10))
        m2.R = fill(0.0, (3,10))
        @test_throws Exception m2.V = fill(0.0, (10,3))
        @test_throws Exception m2.V = fill(nothing, (3,10))
        m2.V = fill(0.0, (3,10))
        @test_throws Exception m2.F = fill(0.0, (10,3))
        @test_throws Exception m2.F = fill(nothing, (3,10))
        m2.F = fill(0.0, (3,10))
        @test_throws Exception m2.ids = fill(0, (10,3))
        @test_throws Exception m2.ids = fill(nothing, 10)
        m2.ids = 1:10
        @test_throws Exception m2.resids = fill(0, (10,3))
        @test_throws Exception m2.resids = fill(nothing, 10)
        m2.resids = fill(1, 10)
        @test_throws Exception m2.names = fill(nothing, 10)
        m2.names = fill("cat", 10)
        @test_throws Exception m2.resnames = fill(nothing, 10)
        m2.resnames = fill("kitten", 10)
        @test_throws Exception m2.chainids = fill(nothing, 10)
        m2.chainids = fill("A", 10)
        @test_throws Exception m2.invalid = fill(nothing, 10)
        @test haskey(m2, :resnames) == true
        @test haskey(m2, :tapir) == false
        @test getkey(m2, :resnames, nothing) == :resnames
        @test getkey(m2, :tapir, nothing) == nothing
        @test get(m2, :resnames, nothing) == fill("kitten", 10)
        @test get(m2, :tapir, nothing) == nothing
        @test_throws Exception m2.tapir = fill("kitten", 10)
        @test get!(m2, :resnames, nothing) == fill("kitten", 10)
        @test get!(m2, :occupancies, fill(1.0, 10)) == fill(1.0, 10)
        get!(m2, :R, undef)
    end

    @testset "Particle" begin
        m1 = MolecularModel(10)
        for (i, p) in enumerate(m1)
            @test p.i == i
        end
        m1.R = fill(0.0, (3,10))
        m1.ids = 1:10
        @test m1[1].R == [0.0, 0.0, 0.0]
        @test m1[2].id == 2
        show(IOBuffer(), m1[3])
        p1 = m1[4]
        @test length(keys(p1)) == 2
        @test_throws Exception p1.name = "CA"
        @test haskey(p1, :id) == true
        @test haskey(p1, :tapir) == false
        @test getkey(p1, :id, nothing) == :id
        @test getkey(p1, :tapir, nothing) == nothing
        @test get(p1, :id, nothing) == 4
        @test get(p1, :tapir, nothing) == nothing
        @test p1.id == 4
        p1.id = 3
    end

    @testset "Views" begin
        m1 = MolecularModel(10)
        v1 = view(m1, 1:5)
        @test parent(v1) == m1
        @test parentindices(v1) == 1:5
        m1.R = fill(0.0, (3,10))
        m1.ids = 1:10
        @test v1.ids == 1:5
        @test v1.R == fill(0.0, (3,5))
        fill!(v1.R, 1.0)
        @test v1.R[:,1:5] == fill(1.0, (3,5))
        @test length(keys(v1)) == 2
        @test_throws Exception v1.resids = 1:5
        @test haskey(v1, :ids) == true
        @test haskey(v1, :tapir) == false
        @test getkey(v1, :ids, nothing) == :ids
        @test getkey(v1, :tapir, nothing) == nothing
        @test get(v1, :ids, nothing) == 1:5
        @test get(v1, :tapir, nothing) == nothing
    end

    @testset "Collection" begin

        @testset "Copy" begin
            m1 = MolecularModel(20)
            m1.ids = 1:20
            m1.header.title = "capybara"
            m2 = MolecularModel(m1)
            @test m2.header.title == "capybara"
            m3 = MolecularModel(view(m1, 6:10))
            @test m3.ids == 6:10
            @test m3.header.title == "capybara"
        end

        @testset "Delete" begin
            m1 = MolecularModel(10)
            m1.ids = 1:10
            @test_throws Exception deleteat!(m, 11)
            m1.topology = []
            pair!(m1.topology, (3, 4))
            pair!(m1.topology, (5, 6))
            deleteat!(m1, 3)
            @test pairs(m1.topology) == [(4,5)]
            @test m1.ids == [1:2..., 4:10...]
            deleteat!(m1, 3:5)
            @test isempty(pairs(m1.topology))
            @test m1.ids == [1:2..., 7:10...]
            pair!(m1.topology, (1, 2))
            empty!(m1)
            @test isempty(m1)
            @test isempty(pairs(m1.topology))
        end

        @testset "Splice" begin
            m1 = MolecularModel(10)
            m1.ids = 1:10
            m1.names = fill("A", 10)
            m2 = MolecularModel(10)
            m2.ids = 11:20
            splice!(m1, 2:3)
            @test length(m1) == 8
            @test m1.ids == [1, 4:10...]
            @test_throws Exception insert!(m1, 2, m2)
            m2.names = fill("B", 10)
            insert!(m1, 2, m2)
            @test length(m1) == 18
            @test m1.ids == [1, 11:20..., 4:10...]
            m1.topology = []
            pair!(m1.topology, (1, 2))
            m2.topology = []
            prepend!(m1, m2[1])
            @test length(m1) == 19
            @test ! ((1, 2) in m1.topology)
            @test (2, 3) in m1.topology
            pair!(m2.topology, (1, 2))
            pair!(m2.topology, (2, 3))
            append!(m1, view(m2, 1:5))
            @test length(m1) == 24
            @test (20, 21) in m1.topology
            @test (21, 22) in m1.topology
        end

        @testset "Resize" begin
            m1 = MolecularModel()
            resize!(m1, 10)
            m1.topology = []
            pair!(m1.topology, (6, 7))
            m1.ids = 1:10
            resize!(m1, 20)
            @test length(m1) == 20
            @test length(pairs(m1.topology)) == 1
            @test m1.ids[5] == 5
            resize!(m1, 5)
            @test length(m1) == 5
            @test isempty(pairs(m1.topology))
        end

    end

    @testset "Group splits" begin

        @testset "Chain" begin
            m1 = MolecularModel(10)
            m1.chainids =
                    [fill("A", 3)..., fill("B", 2)..., fill("C", 4)..., "D"]
            h = hierarchy(m1)
            @test [length(chain) for chain in h.chains] == [3,2,4,1]
            @test parentindices(h.chains[1]) == 1:3
            @test parentindices(h.chains[2]) == 4:5
            @test parentindices(h.chains[3]) == 6:9
            @test parentindices(h.chains[4]) == 10:10
            @test chainindices(m1, 1) == 1:3
            @test chainindices(m1, 2) == 1:3
            @test chainindices(m1, 3) == 1:3
            @test chainindices(m1, 4) == 4:5
            @test chainindices(m1, 5) == 4:5
            @test chainindices(m1, 6) == 6:9
            @test chainindices(m1, 7) == 6:9
            @test chainindices(m1, 8) == 6:9
            @test chainindices(m1, 9) == 6:9
            @test chainindices(m1, 10) == 10:10
        end

        @testset "Residue" begin
            m1 = MolecularModel(20)
            m1.resids = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5]
            m1.resnames = fill("WAT", 20)
            @test [length(res) for res in eachresidue(m1)] == [4,4,4,4,4]
            m1.chainids = [fill("A", 10)..., fill("B", 10)...]
            h = hierarchy(m1)
            @test [length(res) for res in h.residues] == [4,4,2,2,4,4]
            @test parentindices(h.residues[1]) == 1:4
            @test parentindices(h.residues[2]) == 5:8
            @test parentindices(h.residues[3]) == 9:10
            @test parentindices(h.residues[4]) == 11:12
            @test parentindices(h.residues[5]) == 13:16
            @test parentindices(h.residues[6]) == 17:20
            @test residueindices(m1, 1) == 1:4
            @test residueindices(m1, 2) == 1:4
            @test residueindices(m1, 3) == 1:4
            @test residueindices(m1, 4) == 1:4
            @test residueindices(m1, 5) == 5:8
            @test residueindices(m1, 16) == 13:16
            @test residueindices(m1, 17) == 17:20
            @test residueindices(m1, 18) == 17:20
            @test residueindices(m1, 19) == 17:20
            @test residueindices(m1, 20) == 17:20
        end

        @testset "Fragment" begin
            m1 = MolecularModel(10)
            m1.topology = []
            pair!(m1.topology, (1, 2), (2, 3), (4, 5), (4, 6), (4, 7), (4, 8))
            frags = eachfragment(m1)
            @test [length(frag) for frag in frags] == [3,5,1,1]
            @test parentindices(frags[1]) == 1:3
            @test parentindices(frags[2]) == 4:8
            @test parentindices(frags[3]) == 9:9
            @test parentindices(frags[4]) == 10:10
            @test fragmentindices(m1, 1) == 1:3
            @test fragmentindices(m1, 2) == 1:3
            @test fragmentindices(m1, 3) == 1:3
            @test fragmentindices(m1, 4) == 4:8
            @test fragmentindices(m1, 5) == 4:8
            @test fragmentindices(m1, 6) == 4:8
            @test fragmentindices(m1, 7) == 4:8
            @test fragmentindices(m1, 8) == 4:8
            @test fragmentindices(m1, 9) == 9:9
            @test fragmentindices(m1, 10) == 10:10
        end

    end

end # @testset
