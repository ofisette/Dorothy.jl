# Model/chain/residue and model/fragment hierarchies

module Hierarchies

using ..Dorothy
using ..Dorothy.Graphs
using ..Dorothy.Multicollections
using ..Dorothy.Utils

export
		mcrptree, mcrp, chains, residues, mfptree, mfp, fragments,
		chainat, residueat, fragmentat

function mcrptree(model::ParticleCollection)
	n = length(model)
	chainids = get(model, :chainids, Repeated("", n))
	resids = get(model, :resids, Repeated(0, n))
	mcrptree(chainids, resids)
end

function mcrptree(chainids::AbstractVector{<:AbstractString},
		resids::AbstractVector{<:Integer})
	n = length(chainids)
	@boundscheck length(resids) == n ||
			error("size mismatch between chainid and resid arrays")
	model = Vector{Vector{Int}}[]
	if n > 0
		residue = [1]
		chain = [residue]
		push!(model, chain)
		lastchainid = chainids[1]
		lastresid = resids[1]
		for i = 2:n
			thischainid = chainids[i]
			thisresid = resids[i]
			if thischainid != lastchainid
				residue = [i]
				chain = [residue]
				push!(model, chain)
			elseif thisresid != lastresid
				residue = [i]
				push!(chain, residue)
			else
				push!(residue, i)
			end
			lastchainid = thischainid
			lastresid = thisresid
		end
	end
	model
end

mcrp(model::ParticleCollection) = mcrp(mcrptree(model), length(model))

mcrp(tree::Vector{Vector{Vector{Int}}}, n::Integer) = H3Hierarchy(tree, n)

chains(model::ParticleCollection) = chains(model, mcrp(model))

chains(model::ParticleCollection, mcrp::H3Hierarchy) =
		H3IteratorH2(model, mcrp)

residues(model::ParticleCollection) = residues(model, mcrp(model))

residues(model::ParticleCollection, mcrp::H3Hierarchy) =
		H2Iterator(model, mcrp.flath2)

residues(model::ParticleCollection, chainindex::Integer) =
		residues(model, chainindex, mcrp(model))

residues(model::ParticleCollection, ichain::Integer, mcrp::H3Hierarchy) =
		H3IteratorH1(model, ichain, mcrp)

mfptree(model::ParticleCollection) =
		mfptree(get(model, :topology, Graph(length(model))))

function mfptree(G::AbstractGraph)
	n = length(G)
	model = Vector{Int}[]
	if n > 0
		visited = falses(n)
		i = 1
		while true
			fragment = sort!(connected!(Int[], G, i, visited))
			push!(model, fragment)
			i = findnext(!, visited, i+1)
			if i == nothing
				break
			end
		end
	end
	model
end

mfp(model::ParticleCollection) = mfp(mfptree(model), length(model))

mfp(tree::Vector{Vector{Int}}, n::Integer) = H2Hierarchy(tree, n)

fragments(model::ParticleCollection) = fragments(model, mfp(model))

fragments(model::ParticleCollection, mfp::H2Hierarchy) = H2Iterator(model, mfp)

function chainat(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	chainids = get(model, :chainids, Repeated("", length(model)))
	view(model, chainat(chainids, i))
end

function chainat(chainids::AbstractArray{<:AbstractString}, i::Integer)
	@boundscheck checkbounds(chainids, i)
	n = length(chainids)
	thischainid = chainids[i]
	firsti = i - 1
	while firsti > 0 && chainids[firsti] == thischainid
		firsti -= 1
	end
	firsti += 1
	lasti = i + 1
	while lasti <= n && chainids[lasti] == thischainid
		lasti += 1
	end
	lasti -= 1
	firsti:lasti
end

function residueat(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	n = length(model)
	chainids = get(model, :chainids, Repeated("", n))
	resids = get(model, :resids, Repeated(0, n))
	view(model, residueat(chainids, resids, i))
end

function residueat(chainids::AbstractArray{<:AbstractString},
		resids::AbstractArray{<:Integer}, i::Integer)
	n = length(chainids)
	@boundscheck begin
		length(resids) == n ||
				error("size mismatch between chainid and resid arrays")
		checkbounds(chainids, i)
	end
	thischainid = chainids[i]
	thisresid = resids[i]
	firsti = i - 1
	while firsti > 0 && chainids[firsti] == thischainid &&
			resids[firsti] == thisresid
		firsti -= 1
	end
	firsti += 1
	lasti = i + 1
	while lasti <= n && chainids[lasti] == thischainid &&
			resids[lasti] == thisresid
		lasti += 1
	end
	lasti -= 1
	firsti:lasti
end

function fragmentat(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	topology = get(model, :topology, Graph(length(model)))
	view(model, fragmentat(topology, i))
end

fragmentat(G::AbstractGraph, i::Integer) = sort!(connected(G, i))

end # module
