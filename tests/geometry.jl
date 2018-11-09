using Test
using Dorothy
using Dorothy.Utils
using Dorothy.Geometry
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Geometry" begin

    @testset "Distances" begin
        P1 = [0.0, 0.0, 0.0]
        P2 = [1.0, 1.0, 1.0]
        P3 = [2.0, 3.0, 4.0]
        @test sqnorm(P1) == 0.0
        @test sqnorm(P2) == 3.0
        @test sqnorm(P3) == 29.0
        @test dist(P1, P1) == 0.0
        @test dist(P1, P2) == dist(P2, P1) == √3.0
        @test dist(P1, P3) == dist(P3, P1) == √29.0
        @test sqdist(P1, P1) == 0.0
        @test sqdist(P1, P2) == sqdist(P2, P1) == 3.0
        @test sqdist(P1, P3) == sqdist(P3, P1) == 29.0
    end

    @testset "Angles" begin
        P1 = [1.0, 0.0, 0.0]
        P2 = [0.0, 0.0, 0.0]
        P3 = [0.0, 1.0, 0.0]
        P4 = [1.0, 1.0, 1.0]
        @test angle(P1, P2, P1) == 0.0
        @test angle(P1, P2, P3) ≈ angle(P3, P2, P1) ≈ τ/4
        @test angle(P2, P1, P3) ≈ angle(P2, P3, P1) ≈ τ/8
        @test angle(P1, P2, P4) ≈ angle(P4, P2, P1) ≈ 0.9553166181245092
    end

    @testset "Dihedrals" begin
        P1 = [1.0, 0.0, 0.0]
        P2 = [0.0, 0.0, 0.0]
        P3 = [0.0, 0.0, 1.0]
        P4 = [1.0, 0.0, 1.0]
        P5 = [-1.0, 0.0, 1.0]
        P6 = [0.0, 1.0, 1.0]
        P7 = [0.0, -1.0, 1.0]
        @test dihedral(P1, P2, P3, P4) == dihedral(P4, P3, P2, P1) == 0.0
        @test dihedral(P1, P2, P3, P5) == dihedral(P5, P3, P2, P1) == τ/2
        @test dihedral(P1, P2, P3, P6) == dihedral(P6, P3, P2, P1) == τ/4
        @test dihedral(P1, P2, P3, P7) == dihedral(P7, P3, P2, P1) == -τ/4
    end

    @testset "COM" begin
        P1 = [1.0, 0.0, 0.0]
        P2 = [0.0, 1.0, 0.0]
        P3 = [0.0, 0.0, 1.0]
        R = [P1 P2 P3]
        @test cog(R) == [1/3, 1/3, 1/3]
    end

    @testset "Transformations" begin
        P1 = [1.0, 0.0, 0.0]
        @test translation(x=1.0) * P1 == [2.0, 0.0, 0.0]
        rotatez4 = rotation(τ/4, :z)
        @test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
        rotatez4 = rotation(τ/4, [0.0, 0.0, 1.0])
        @test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
        rotatez4 = rotation(τ/4, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        @test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
        rotatez4 = rotation(τ/4, [0.0, 0.0, 10.0] - [0.0, 0.0, 1.0])
        @test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
        scale2 = scaling(2.0)
        @test scale2 * P1 == [2.0, 0.0, 0.0]
        @test translation(1.0, 0.0, 1.0) * scaling(2.0) * rotation(τ/4, :z) *
                P1 ≈ [1.0, 2.0, 1.0]
    end

    @testset "Superposition" begin
        m1 = readf("$(datapath)/1BTL.pdb")
        m2 = MolecularModel(m1)
        @test rmsd(m1.R, m2.R) ≈ 0.0
        m1.R = translation(25.0, 30.0, -1.0) * rotation(τ/4, :z) *
                translation(-5.0, 3.0, 1.0) * m2.R
        @test rmsd(m1.R, m2.R) ≈ 36.167398917591804
        tr = superposition(m1.R, m2.R)
        m1.R .= tr * m1.R
        @test rmsd(m1.R, m2.R) + 1.0 ≈ 1.0
    end

    @testset "Fit line" begin
        m = readf("$(datapath)/1BTL.pdb")
        h1 = view(m, 1:120) # First α helix
        V, comR = fitline(h1.R)
        @test V ≈ [0.7767115502007541, 0.3946883305177638, 0.4908566894092842]
        @test comR ≈ [12.60863333333333, 7.342849999999999, 12.188675000000002]
    end

    @testset "Fit plane" begin
        m = readf("$(datapath)/1BTL.pdb")
        h18 = view(m, [1:120; 1866:2032]) # First and last α helices
        V, comR = fitplane(h18.R)
        @test V ≈ [-0.02187509375361663, 0.815568953125186, -0.578246282280793]
        @test comR ≈ [9.50672822299652, 10.308843205574918, 15.995000000000005]
    end

    @testset "PBC" begin
        c1 = pbccell(30.0)
        @test c1[1,1] == 30.0
        @test c1[1,2] + 1 ≈ 1.0
        @test c1 == pbccell(30.0, 30.0)
        @test c1 == pbccell(30.0, 30.0, 30.0)
        @test c1 != pbccell(30.0, 25.0)
        c2 = pbccell(pbcbox(c1)...)
        @test c2 == c1
        @test iscubic(pbcbox(c1)...)
        @test volume(c1) == 30.0^3
        @test volume(pbcbox(c1)...) == 30.0^3
    end

    @testset "Wrap" begin
        m1 = readf("$(datapath)/1BTL.pdb")
        @test count(x -> (x < 0.0), m1.R) == 1581 # Coords outside the cell.
        wrappos!(m1.R, m1.header.cell)
        @test count(x -> (x < 0.0), m1.R) == 0
        m2 = readf("$(datapath)/1BTL.pdb")
        m2.R[1,:] .-= 6.0
        m2.R[2,:] .-= 3.0
        cogR = cog(m2.R)
        @test cogR[1] < 0.0
        @test cogR[2] < 0.0
        @test cogR[3] > 0.0
        wrapcog!(m2.R, m2.header.cell)
        cogR = cog(m2.R)
        @test cogR[1] > 0.0
        @test cogR[2] > 0.0
        @test cogR[3] > 0.0
    end

    @testset "Unwrap" begin
        m1 = readf("$(datapath)/1BTL.pdb")
        m1.header.cell = pbccell(48.0, 64.0, 91.0)
        protein = view(m1, 1:2032)
        prewrapdim = dim(extent(protein.R))
        wrappos!(m1.R, m1.header.cell)
        @test all(dim(extent(protein.R)) .>= prewrapdim)
        Dkmax = unwrappedmaxdim(protein.header.cell)
        unwrapbydim!(protein.R, protein.header.cell, Dkmax)
        unwrappedext = extent(protein.R)
        @test dim(unwrappedext) ≈ prewrapdim
        @test all(minimum(unwrappedext) .>= 0.0)
    end

    #=
    @testset "Compact" begin
        # Get a small rhdo box and compactify it. Do the same to a full molecule
        # such as a protein.
        m1 = readf("$(datapath)/MHC.pdb")
        compactpos!(m1.R, m1.header.cell)
        writef("$(datapath)/test.pdb", m1)
    end
    =#

end # @testset
