using Test
using Dorothy

@DorothyAll()

@testset "Topology" begin

	@testset "Inference" begin
		m1 = readf("$(Dorothy.datapath)/1BTL.gro")
		topo = infertopology!(m1)
		@test topo[1] == [2]
		@test topo[2] == [1,3,5]
		@test topo[3] == [2,4,11]
		@test topo[4] == [3]
		m2 = readf("$(Dorothy.datapath)/MHC.pdb")
		topo = infertopology!(m2)
		@test topo[1] == [2,3,4,5,6,7]
		@test topo[2] == [1,3,5,6,7]
		@test topo[3] == [1,2,4,5,6,7]
		@test topo[4] == [1,3]
		@test topo[6540] == [6541, 6542]
		@test length(fragments(view(m2, Protein))) == 3
	end

end # @testset
