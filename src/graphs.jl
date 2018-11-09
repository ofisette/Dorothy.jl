# Simple graph types

module Graphs

using .Iterators: repeated
using ..Dorothy.Utils

export
		AbstractGraph, Graph, FixedGraph, neighbors!, pairs, connected,
		connected!, isisolated, pair!, unpair!, isolate!

abstract type AbstractGraph <: AbstractVector{Vector{Int}} end

struct Graph <: AbstractGraph
	neighbors::Vector{Vector{Int}}

	function Graph(n::Integer)
		n >= 0 || error("expected positive number of indices")
		new(fill(emptyindices, n))
	end
end

struct FixedGraph <: AbstractGraph
	G::Graph
end

struct GraphView{T<:AbstractVector{<:Integer}} <: AbstractGraph
	G::Graph
	I::T

	@inline function GraphView{T}(G::Graph, I::T) where
			{T<:AbstractVector{<:Integer}}
		@boundscheck checkindexseries(G, I)
		new(G, I)
	end
end

function Graph(G::AbstractGraph)
	G2 = Graph(length(G))
	for (i, j) in pairs(G)
		pair!(G2, (i, j))
	end
	G2
end

@inline GraphView(G::Graph, I::T) where {T<:AbstractVector{<:Integer}} =
		GraphView{T}(G, I)

@inline GraphView(G::FixedGraph, I::T) where {T<:AbstractVector{<:Integer}} =
		GraphView{T}(G.G, I)

@inline function GraphView(G::GraphView, I::AbstractVector{<:Integer})
	@boundscheck checkindexseries(G, I)
	@inbounds GraphView(G.G, G.I[I])
end

Base.size(G::Graph) = (length(G.neighbors),)

Base.size(G::FixedGraph) = (length(G.G),)

Base.size(G::GraphView) = (length(G.I),)

Base.getindex(G::AbstractGraph, i::Int) = neighbors!(Int[], G, i)

Base.getindex(G::Graph, i::Int) = copy(G.neighbors[i])

Base.IndexStyle(::Type{<:AbstractGraph}) = IndexLinear()

Base.parent(G::Graph) = G

Base.parent(G::FixedGraph) = parent(G.G)

Base.parent(G::GraphView) = G.G

Base.parentindices(G::Graph) = 1:length(G)

Base.parentindices(G::FixedGraph) = parentindices(G.G)

Base.parentindices(G::GraphView) = G.I

Base.view(G::AbstractGraph, I::AbstractVector{<:Integer}) = GraphView(G, I)

Base.view(G::AbstractGraph, i::Integer) = GraphView(G, [i])

Base.view(G::Graph, I::AbstractVector{<:Integer}) = GraphView(G, copy(I))

@inline function Base.in((i, j)::Tuple{Integer,Integer}, G::Graph)
	@boundscheck begin
		checkbounds(G, i)
		checkbounds(G, j)
	end
	j in G.neighbors[i]
end

@inline Base.in((i, j)::Tuple{Integer,Integer}, G::FixedGraph) = in((i, j), G.G)

@inline function Base.in((i, j)::Tuple{Integer,Integer}, G::GraphView)
	@boundscheck begin
		checkbounds(G, i)
		checkbounds(G, j)
	end
	(G.I[i], G.I[j]) in G.G
end

@inline function neighbors!(dest::AbstractArray{<:Integer}, G::Graph,
		i::Integer)
	@boundscheck checkbounds(G, i)
	empty!(dest)
	append!(dest, G.neighbors[i])
end

@inline neighbors!(dest::AbstractArray{<:Integer}, G::FixedGraph, i::Integer) =
		neighbors!(dest, G.G, i)

@inline function neighbors!(dest::AbstractArray{<:Integer}, G::GraphView,
		i::Integer)
	@boundscheck checkbounds(G, i)
	neighbors!(dest, G.G, G.I[i])
	for j in reverse(eachindex(dest))
		r = searchsorted(G.I, dest[j])
		if isempty(r)
			deleteat!(dest, j)
		else
			dest[j] = r[1]
		end
	end
	dest
end

function Base.pairs(G::Graph)
	pairs = sizehint!(Tuple{Int,Int}[], length(G))
	for i in eachindex(G.neighbors)
		for j in G.neighbors[i]
			if j > i
				push!(pairs, (i, j))
			end
		end
	end
	pairs
end

Base.pairs(G::FixedGraph) = pairs(G.G)

function Base.pairs(G::GraphView)
	pairs = sizehint!(Tuple{Int,Int}[], length(G))
	tmp = Int[]
	for i in 1:length(G)
		for j in @inbounds neighbors!(tmp, G, i)
			if j > i
				push!(pairs, (i, j))
			end
		end
	end
	pairs
end

connected(G::AbstractGraph, i::Integer) = connected!(Int[], G, i)

function connected!(dest::AbstractVector{<:Integer}, G::AbstractGraph,
		i::Integer, skip::AbstractVector{<:Bool} = falses(length(G)))
	@boundscheck checkbounds(G, i)
	empty!(dest)
	tmp = Int[]
	stack = [i]
	while ! isempty(stack)
		j = pop!(stack)
		if ! skip[j]
			push!(dest, j)
			append!(stack, neighbors!(tmp, G, j))
			skip[j] = true
		end
	end
	dest
end

function isisolated(G::AbstractGraph, I::AbstractArray{<:Integer})
	@boundscheck checkbounds(G, I)
	J = Int[]
	for i in I
		for j in neighbors!(J, G, i)
			if ! (j in I)
				return false
			end
		end
	end
	true
end

@inline function isisolated(G::Graph, i::Integer)
	@boundscheck checkbounds(G, i)
	isempty(G.neighbors[i])
end

@inline isisolated(G::FixedGraph, i::Integer) = isisolated(G.G, i)

@inline function isisolated(G::GraphView, i::Integer)
	@boundscheck checkbounds(G, i)
	@inbounds isisolated(G.G, G.I[i])
end

function pair!(G::Graph, pairs::Tuple{Integer,Integer}...)
	@boundscheck begin
		for (i, j) in pairs
			checkbounds(G, i)
			checkbounds(G, j)
			i != j || error("expected non-identical indices")
		end
	end
	for (i, j) in pairs
		J = G.neighbors[i]
		r = searchsorted(J, j)
		if isempty(r)
			if isempty(J)
				G.neighbors[i] = [j]
			else
				insert!(J, first(r), j)
			end
			I = G.neighbors[j]
			if isempty(I)
				G.neighbors[j] = [i]
			else
				insert!(I, searchsortedfirst(I, i), i)
			end
		end
	end
	G
end

pair!(G::FixedGraph, pairs::Tuple{Integer,Integer}...) = pair!(G.G, pairs...)

function pair!(G::GraphView, pairs::Tuple{Integer,Integer}...)
	@boundscheck begin
		for (i, j) in pairs
			checkbounds(G, i)
			checkbounds(G, j)
			i != j || error("expected non-identical indices")
		end
	end
	@inbounds pair!(G.G, [(G.I[i], G.I[j]) for (i, j) in pairs]...)
end

function unpair!(G::Graph, pairs::Tuple{Integer,Integer}...)
	@boundscheck begin
		for (i, j) in pairs
			checkbounds(G, i)
			checkbounds(G, j)
			i != j || error("expected non-identical indices")
		end
	end
	for (i, j) in pairs
		J = G.neighbors[i]
		r = searchsorted(J, j)
		if ! isempty(r)
			if length(J) == 1
				G.neighbors[i] = emptyindices
			else
				deleteat!(J, first(r))
			end
			I = G.neighbors[j]
			if length(I) == 1
				G.neighbors[j] = emptyindices
			else
				deleteat!(I, searchsortedfirst(I, i))
			end
		end
	end
	G
end

unpair!(G::FixedGraph, pairs::Tuple{Integer,Integer}...) =
		unpair!(G.G, pairs...)

function unpair!(G::GraphView, pairs::Tuple{Integer,Integer}...)
	@boundscheck begin
		for (i, j) in pairs
			checkbounds(G, i)
			checkbounds(G, j)
			i != j || error("expected non-identical indices")
		end
	end
	@inbounds unpair!(G.G, [(G.I[i], G.I[j]) for (i, j) in pairs]...)
end

function isolate!(G::Graph, i::Integer)
	@boundscheck checkbounds(G, i)
	for j in G.neighbors[i]
		I = G.neighbors[j]
		if length(I) == 1
			G.neighbors[j] = emptyindices
		else
			deleteat!(I, searchsortedfirst(I, i))
		end
	end
	G.neighbors[i] = emptyindices
	G
end

@inline isolate!(G::FixedGraph, i::Integer) = isolate!(G.G, i)

@inline function isolate!(G::GraphView, i::Integer)
	@boundscheck checkbounds(G, i)
	for j in @inbounds G[i]
		unpair!(G, (i, j))
	end
	G
end

function isolate!(G::AbstractGraph, I::AbstractVector{<:Integer})
	@boundscheck checkbounds(G, I)
	J = Int[]
	for i in I
		for j in neighbors!(J, G, i)
			if ! (j in I)
				unpair!(G, (i, j))
			end
		end
	end
	G
end

function Base.:(==)(G1::AbstractGraph, G2::AbstractGraph)
	if length(G1) == length(G2)
		J1 = Int[]
		J2 = Int[]
		for i in eachindex(G1)
			J1 = neighbors!(G1, i)
			J2 = neighbors!(G2, i)
			if J1 != J2
				return false
			end
		end
		true
	else
		false
	end
end

function Base.merge!(G::AbstractGraph, src::AbstractGraph)
	@boundscheck length(src) <= length(G) ||
			error("cannot merge source graph larger than destination")
	pair!(G, pairs(src)...)
	G
end

function Base.resize!(G::Graph, n::Integer)
	@boundscheck n >= 0 || error("expected positive length")
	nprev = length(G)
	dn = n - nprev
	if dn > 0
		resize!(G.neighbors, n)
		for i = nprev+1:n
			G.neighbors[i] = emptyindices
		end
	elseif dn < 0
		for i = n+1:nprev
			for j in G.neighbors[i]
				I = G.neighbors[j]
				if length(I) == 1
					G.neighbors[j] = emptyindices
				else
					deleteat!(I, searchsortedfirst(I, i))
				end
			end
			G.neighbors[i] = emptyindices
		end
		resize!(G.neighbors, n)
	end
	G
end

function Base.empty!(G::Graph)
	empty!(G.neighbors)
	G
end

function remap!(G::Graph, map::AbstractVector{<:Integer})
	newneighbors = Int[]
	for (i0, i1) in enumerate(map)
		if i1 != 0
			empty!(newneighbors)
			for j0 in G.neighbors[i0]
				j1 = map[j0]
				if j1 != 0
					push!(newneighbors, j1)
				end
			end
			if isempty(newneighbors)
				G.neighbors[i0] = emptyindices
			else
				G.neighbors[i0] = copy(newneighbors)
			end
		end
	end
end

Base.deleteat!(G::Graph, i::Integer) = deleteat!(G, [i])

@inline function Base.deleteat!(G::Graph, I::AbstractArray{<:Integer})
	@boundscheck checkindexseries(G, I)
	remap!(G, deleteat_map(length(G), I))
	deleteat!(G.neighbors, I)
	G
end

Base.splice!(G::Graph, i::Integer) = deleteat!(G, i)

Base.splice!(G::Graph, i::Integer, replacement::AbstractGraph) =
		splice!(G, i:i, replacement)

Base.splice!(G::Graph, range::UnitRange{<:Integer}) = deleteat!(G, range)

function Base.splice!(G::Graph, range::UnitRange{<:Integer},
		replacement::AbstractGraph)
	@boundscheck begin
		checkbounds(G, range)
		parent(replacement) !== G ||
				error("replacement must belong to a different graph")
	end
	spliced = Graph(view(G, range))
	nreplacement = length(replacement)
	if nreplacement == 0
		@inbounds deleteat!(G, range)
	else
		remap!(G, splice_map(length(G), range, nreplacement))
		splice!(G.neighbors, range, repeated(emptyindices, nreplacement))
		offset = first(range) - 1
		for (i, j) in pairs(replacement)
			pair!(G, (i+offset, j+offset))
		end
	end
	spliced
end

function Base.append!(G::Graph, src::AbstractGraph)
	n = length(G)
	splice!(G, n+1:n, src)
	G
end

function Base.prepend!(G::Graph, src::AbstractGraph)
	splice!(G, 1:0, src)
	G
end

function Base.insert!(G::Graph, i::Integer, src::AbstractGraph)
	splice!(G, i:i-1, src)
	G
end

end # module

#= Crusty old stuff from before this Graphs module

struct PathAtBondDistIterator{T<:ParticleCollection}
	model::T
	start::Int
	dist::Int

	function PathAtBondDistIterator{T}(model::T, start::Integer,
			dist::Integer) where {T<:ParticleCollection}
		@boundscheck begin
			checkbounds(model, start)
			dist > 0 || error("expected strictly positive bond distance")
		end
		new(model, start, dist)
	end
end

Base.iterate(iter::PathAtBondDistIterator) = iterate(iter, [[iter.start]])

function Base.iterate(iter::PathAtBondDistIterator, stack::Vector{Vector{Int}})
	tmp = Int[]
	while length(stack) != iter.dist + 1 || isempty(stack[end])
		while isempty(stack[end])
			pop!(stack)
			if isempty(stack)
				return
			else
				pop!(stack[end])
			end
		end
		newlevel = Int[]
		for i in bondedpartners!(tmp, iter.model, stack[end][end])
			visited = false
			for level in stack
				if i == level[end]
					visited = true
					break
				end
			end
			if ! visited
				push!(newlevel, i)
			end
		end
		push!(stack, newlevel)
	end
	path = [level[end] for level in stack]
	pop!(stack[end])
	path, stack
end

Base.IteratorSize(::PathAtBondDistIterator) = Base.SizeUnknown()

Base.eltype(::PathAtBondDistIterator) = Vector{Int}

distbondedpaths(model::T, start::Integer, dist::Integer) where
		{T<:ParticleCollection} = PathAtBondDistIterator{T}(model, start, dist)

function bondedpaths(model::ParticleCollection, i::Integer, j::Integer,
		maxdist::Integer = length(model) - 1)
	@boundscheck begin
		checkbounds(model, i)
		checkbounds(model, j)
		i != j || error("cannot find paths between a particle and itself")
		maxdist > 0 || error("expected strictly positive maximum bond distance")
		maxdist < length(model) || error("maximum bond distance must be" *
				"smaller than number of particles")
	end
	paths = Vector{Int}[]
	for dist = 1:maxdist
		for path in distbondedpaths(model, i, dist)
			if path[end] == j
				push!(paths, path)
			end
		end
	end
	paths
end

function bondedpath(model::ParticleCollection, i::Integer, j::Integer,
		maxdist::Integer = length(model) - 1)
	@boundscheck begin
		checkbounds(model, i)
		checkbounds(model, j)
		i != j || error("cannot find path between a particle and itself")
		maxdist > 0 || error("expected strictly positive maximum bond distance")
		maxdist < length(model) || error("maximum bond distance must be" *
				"smaller than number of particles")
	end
	for dist = 1:maxdist
		for path in distbondedpaths(model, i, dist)
			if path[end] == j
				return path
			end
		end
	end
	nothing
end
=#
