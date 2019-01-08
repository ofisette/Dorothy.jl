using Test
using Dorothy
using Dorothy.Utils
using Dorothy.Geometry
using Formats
using FormatStreams

using Dorothy.Selectors

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
		R = Vector3D[P1, P2, P3]
		@test cog(R) == [1/3, 1/3, 1/3]
	end

	@testset "Transformations" begin
		P1 = [1.0, 0.0, 0.0]
		@test translation(x = 1.0) * P1 == [2.0, 0.0, 0.0]
		rotatez4 = rotation(z = τ/4)
		@test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
		rotatez4 = rotation(τ/4, [0.0, 0.0, 1.0])
		@test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
		rotatez4 = rotation(τ/4, [0.0, 0.0, 1.0], O = [0.0, 0.0, 0.0])
		@test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
		rotatez4 = rotation(τ/4, [0.0, 0.0, 10.0] - [0.0, 0.0, 1.0])
		@test rotatez4 * P1 ≈ [0.0, 1.0, 0.0]
		scale2 = scaling(2.0)
		@test scale2 * P1 == [2.0, 0.0, 0.0]
		@test translation([1.0, 0.0, 1.0]) * scaling(2.0) * rotation(z = τ/4) *
				P1 ≈ [1.0, 2.0, 1.0]
	end

	@testset "Superposition" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		m2 = MolecularModel(m1)
		@test rmsd(m1.R, m2.R) ≈ 0.0
		m1.R .= translation([25.0, 30.0, -1.0]) * rotation(z = τ/4) *
				translation([-5.0, 3.0, 1.0]) .* m2.R
		@test rmsd(m1.R, m2.R) ≈ 36.167398917591804
		T = superposition(m1.R, m2.R)
		m1.R .= T .* m1.R
		@test rmsd(m1.R, m2.R) + 1.0 ≈ 1.0
	end

	@testset "Fit line" begin
		m = readf("$(datapath)/1BTL.pdb")
		guessmasses!(m)
		h1 = view(m, 1:120) # First α helix
		V = fitline(h1.R)
		@test V ≈ [0.776711550200754, 0.3946883305177638, 0.4908566894092842]
		V = fitline(h1.R, h1.masses)
		@test V ≈ [0.783371080605647, 0.3869484031172630, 0.4864161627616218]
	end

	@testset "Fit plane" begin
		m = readf("$(datapath)/1BTL.pdb")
		guessmasses!(m)
		h18 = view(m, [1:120; 1866:2032]) # First and last α helices
		V = fitplane(h18.R)
		@test V ≈ [-0.02187509375361663, 0.815568953125186, -0.578246282280793]
		V = fitplane(h18.R, h18.masses)
		@test V ≈ [-0.02428341812381785, 0.828435699705839, -0.559557510053364]
	end

end # @testset
