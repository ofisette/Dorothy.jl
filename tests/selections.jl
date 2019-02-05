using Test
using Dorothy

@DorothyAll()

datapath = joinpath(@__DIR__, "..", "data")

@testset "Selection" begin

	@testset "1BTL" begin # Small orthorhombic system
		model = readf("$(datapath)/1BTL.pdb")
		sel1 = view(model, Water)
		@test sel1.ids == 2039:2237
		@test length(view(model, Protein)) == 2032
		@test length(view(model, Water)) == 199
		@test length(view(model, Water | Protein)) == 2032 + 199
		@test length(view(model, Water & Protein)) == 0
		@test length(view(model, Index(1, by=Residue))) == 10
		@test length(view(model, Expand(Index(1), by=Residue))) == 10
		@test length(view(model, Index(1, by=Chain))) == 2236
		@test length(view(model, Expand(Index(1), by=Chain))) == 2236
		sel2 = view(model, Water & Within(5.0, of=Nter))
		@test length(sel2) == 7
		sel3 = view(model, Water & Within(45.0, of=Nter))
		@test length(sel3) == 199
	end

	@testset "1BTL-np" begin # Non-periodic system
		model = readf("$(datapath)/1BTL.pdb")
		delete!(model.header, :cell)
		sel1 = view(model, Water & Within(45.0, of=Nter))
		@test length(sel1) == 188
	end

	@testset "MHC" begin # Rhombododecahedral system
		model = readf("$(datapath)/MHC.pdb")
		sel1 = view(model, Water & Within(5.0, of=Protein))
		@test length(sel1) == 6475
		sel2 = view(model, Cter)
		@test length(sel2) == 56
		@test length(residues(sel2)) == 3
		sel3 = view(model, Nter)
		@test length(sel3) == 57
		@test length(residues(sel3)) == 3
		bitvec = falses(length(model))
		for i = 1:6539
			bitvec[i] = true
		end
		sel4 = view(model, Map(bitvec))
		@test length(sel4) == 6539
		sel5 = view(model, Protein)
		@test sel5 == sel4
		@test length(sel5) == 6539
		sel6 = view(model, Index(1:11))
		@test length(sel6) == 11
		@test all(sel6.resnames .== "GLY")
		sel7 = view(model, Index(1))
		@test sel7[].name == "MN1"
		sel8 = view(model, Index(277, by=Residue))
		@test length(sel8) == 27
		@test all(sel8.resnames .== "ILE")
		sel9 = view(model, Index(3, by=Chain))
		@test length(sel9) == 149
		@test all(sel9.chainids .== "C")
		sel10 = view(model, Index(first, ineach=Chain))
		@test length(sel10) == 6
		Dorothy.Topology.infertopology!(model)
		@test length(fragments(view(model, 1:6539))) == 3
		sel11 = view(model, Index(3, by=Fragment))
		@test sel11 == sel9
		@test length(fragments(view(model, Protein))) == 3
		sel12 = view(model, Expand(Index(6539), by=Chain))
		@test sel12 == sel11
		sel13 = view(model, Expand(Index(4657), by=Residue))
		@test sel13 == sel8
		sel14 = view(model, Restrict(Index(6000:7000), by=Chain))
		@test sel14 == sel11
		sel15 = view(model, Restrict(Index(4650:4690), by=Residue))
		@test sel15 == sel8
		sel16 = view(model, Restrict(Index(6541:6565), by=Fragment))
		@test length(sel16) == 3*7
		sel17 = view(model, Heavy & Water)
		@test length(sel17) == 20579
		sel18 = view(model, AcidResidue & Cα)
		sel19 = view(model, BasicResidue & Cα)
		sel20 = view(model, ChargedResidue & Cα)
		sel21 = view(model, PolarResidue & Cα)
		sel22 = view(model, HydrophobicResidue & Cα)
		@test sum(length.([sel18, sel19])) == length(sel20)
		@test sum(length.([sel20, sel21, sel22])) == 276 + 99 + 9
		sel23 = view(model, Backbone)
		@test length(sel23) == 3 * (276 + 99 + 9)
		sel24 = view(model, MainChain)
		@test length(sel24) == 5 * (276 + 99 + 9) + 2*3 + 1*3 - (16 + 5)
		sel25 = view(model, SideChain & MainChain)
		@test length(sel25) == 0
		sel26 = view(model, SideChain | MainChain)
		@test length(sel26) == length(sel5)
	end

	@testset "PLC" begin # Very large orthorhombic system
		model = readf("$(datapath)/PLC.pdb.bz2")
		sel1 = view(model, Water & Expand(Within(3.0, of=Lipid), by=Residue))
		@test length(sel1) == 30474
		sel2 = view(model, Water & Restrict(Within(3.0, of=Lipid), by=Residue))
		@test length(sel2) == 14271
		sel3 = view(model, Name("CA"))
		@test length(sel3) == 4854
		sel4 = view(model, Id(1:10))
		@test length(sel4) == 10
		sel5 = view(model, ResName("LEU"))
		@test length(sel5) == 11708
		sel6 = view(model, ResName("LEU") & Name("CA"))
		@test length(sel6) == 509
		sel7 = view(model, ResName("LEU") & ResId(732))
		@test length(sel7) == 24
		cache = SelectionCache()
		sel8 = view(model, ChainId("X") & Element("Cl"), cache)
		@test length(sel8) == 679
		sel9 = view(model, ChainId("X") & Element("Na"), cache)
		@test length(sel9) == 813
		sel10 = view(model, ChainId("X") & !Element("Na"), cache)
		@test sel10 == sel8
		sel11 = view(model, ChainId("X") & Element("Na") & Element("Cl"))
		@test length(sel11) == 0
		sel12 = view(model, ChainId("X") & (Element("Na") | Element("Cl")))
		@test length(sel12) == 679 + 813
		sel13 = view(model, ChainId("X") & Element(["Na", "Cl"]))
		@test sel13 == sel12
		sel14 = view(model, Mass(mi -> (mi < 2.0)))
		sel15 = view(model, Hydrogen | VSite)
		@test sel14 == sel15
		sel16 = view(model, Occupancy(1.0))
		@test length(sel16) == length(model)
		sel17 = view(model, R([120.340, 167.910, 179.540]))
		@test length(sel17) == 1
		sel18 = view(model, R(Ri -> (Ri.z > 100.0)))
		@test length(sel18) == 1104250
		# Stride test requires binary dependency
		# sel19 = view(model, SS("H"))
		sel20 = view(model, Ion)
		@test sel20 == sel12
		sel21 = view(model, Lipid & Element("P"))
		@test length(sel21) == 1943
	end

end # @testset
