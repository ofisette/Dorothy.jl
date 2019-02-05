using Test
using Dorothy

@DorothyAll()

datapath = joinpath(@__DIR__, "..", "data")

@testset "Graphs" begin

	@testset "Basics" begin
		G = Graph(10)
		@test !((6, 7) in G)
		pair!(G, (6, 7))
		@test (6, 7) in G
		pair!(G, (6, 7))
		@test G[6] == [7]
		pair!(G, (6, 8))
		@test G[6] == [7, 8]
		@test_throws Exception G[11]
		@test_throws Exception (10, 11) in G
		@test_throws Exception pair!(G, (1, 11))
		v1 = view(G, 6:10)
		@test (1, 2) in v1
		@test (1, 3) in v1
		@test !((1, 4) in v1)
		@test v1[1] == [2, 3]
		pair!(v1, (4, 5))
		@test (9, 10) in G
		pair!(v1, (4, 5))
		@test_throws Exception v1[6]
		@test_throws Exception (5, 6) in v1
		@test_throws Exception pair!(v1, (1, 6))
		unpair!(G, (6, 10))
		@test_throws Exception unpair!(G, (9, 11))
		unpair!(v1, (1, 5))
		@test_throws Exception unpair!(v1, (9, 11))
		unpair!(v1, (4, 5))
		unpair!(G, (6, 7))
		@test pairs(G) == [(6, 8)]
		@test pairs(v1) == [(1, 3)]
		pair!(G, (1, 2))
		pair!(G, (1, 3))
		isolate!(G, 1)
		@test isempty(G[1])
		isolate!(v1, 1)
		@test isempty(pairs(v1))
		pair!(G, (1, 4))
		pair!(G, (1, 5))
		for i in 1:length(G)
			isolate!(G, i)
		end
		@test isempty(pairs(G))
	end

	@testset "Connected" begin
		G = Graph(10)
		@test isisolated(G, 1:10)
		pair!(G, (1, 2))
		@test isisolated(G, 1:10)
		@test !isisolated(G, 2:10)
		isolate!(G, 2:10)
		@test isisolated(G, 2:10)
		pair!(G, (1, 2))
		pair!(G, (2, 3))
		pair!(G, (4, 5))
		pair!(G, (4, 6))
		pair!(G, (4, 7))
		pair!(G, (4, 8))
		@test sort!(connected(G, 1)) == [1, 2, 3]
		@test sort!(connected(G, 4)) == [4, 5, 6, 7, 8]
	end

end # @testset
