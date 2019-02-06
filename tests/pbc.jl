using Test
using Dorothy

@DorothyAll()

datapath = joinpath(@__DIR__, "..", "data")

@testset "PBC" begin

	@testset "Cell" begin
		c1 = TriclinicCell(OrthorhombicCell(30.0))
		(a,b,c), (α,β,γ) = pbcgeometry(c1)
		@test [a,b,c] == [30.0,30.0,30.0]
		@test [α,β,γ] == [τ/4,τ/4,τ/4]
		M = pbcmatrix(c1)
		@test M[1,1] == 30.0
		@test M[1,2] + 1 ≈ 1.0
		@test c1 == OrthorhombicCell(30.0, 30.0, 30.0)
		@test c1 == OrthorhombicCell(30.0, 30)
		@test c1 != OrthorhombicCell(30.0, 25.0)
		c2 = OrthorhombicCell(30.0)
		@test c2 == c1
		@test ispbc(c1, cubic)
		@test ispbc(c2, cubic)
		@test volume(c1) == 30.0^3
		@test pbcvolume(pbcgeometry(c1)...) == 30.0^3
	end

	@testset "Wrap" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		@test count(R -> !isinside(R, m1.header.cell), m1.R) == 1268
		m1.R .= wrappos.(m1.R, m1.header.cell)
		@test count(R -> !isinside(R, m1.header.cell), m1.R) == 0
		m2 = readf("$(datapath)/1BTL.pdb")
		translation([-6.0, -3.0, 0.0])(m2.R)
		cogR = cog(m2.R)
		@test cogR.x < 0.0 && cogR.y < 0.0 && cogR.z > 0.0
		wrappos!(m2.R, m2.header.cell)
		cogR = cog(m2.R)
		@test cogR.x > 0.0 && cogR.y > 0.0 && cogR.z > 0.0
	end

	@testset "Unwrap" begin
		m1 = readf("$(datapath)/1BTL.pdb")
		m1.header.cell = OrthorhombicCell([48.0, 64.0, 91.0])
		protein = view(m1, 1:2032)
		cell = protein.header.cell
		prewrapdims = dims(extent(protein.R))
		m1.R .= wrappos.(m1.R, m1.header.cell)
		@test all(dims(extent(protein.R)) .>= prewrapdims)
		unwrap!(protein.R, cell, UnwrapByGap(cell))
		unwrappedext = extent(protein.R)
		@test dims(unwrappedext) ≈ prewrapdims
		@test all(minimum(unwrappedext) .>= 0.0)
	end

	@testset "Compact positions" begin
		model = readf("$(datapath)/MHC.pdb")
		R = model.R
		cell = model.header.cell
		R .= wrappos.(R, cell)
		unwrap!(view(model, Protein).R, cell, UnwrapByGap(cell))
		for residue in residues(view(model, Water))
			unwrap!(residue.R, cell, UnwrapByExtent())
		end
		C = center(cell)
		R .= nearestpos.(R, Ref(C), cell)
		@test round.(dims(extent(R))) ≈ [99.0, 99.0, 137.0]
	end

	@testset "Compact COG" begin
		model = readf("$(datapath)/MHC.pdb")
		R = model.R
		cell = model.header.cell
		C = center(cell)
		nearestpos!(view(model, Protein).R, C, cell)
		for residue in residues(view(model, Water))
			nearestpos!(residue.R, C, cell)
		end
		for atom in view(model, Ion)
			atom.R = nearestpos(atom.R, C, cell)
		end
		@test round.(dims(extent(R))) ≈ [100.0, 100.0, 138.0]
	end

end # @testset
