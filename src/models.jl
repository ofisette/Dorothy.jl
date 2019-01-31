# Types for molecular models, their particles and trajectories

struct MolecularModelHeader <: Header
	D::Dict{Symbol,Any}

	MolecularModelHeader() = new(Dict{Symbol,Any}())
end

function MolecularModelHeader(header::MolecularModelHeader)
	header2 = MolecularModelHeader()
	merge!(header2, header)
end

Headers.headerval(::MolecularModelHeader, ::Val{:title}, v::AbstractString) =
		String(v)

Headers.headerval(::MolecularModelHeader, ::Val{:step}, v::Integer) = Int(v)

Headers.headerval(::MolecularModelHeader, ::Val{:time}, v::Real) = Float64(v)

Headers.headerval(::MolecularModelHeader, ::Val{:lambda}, v::Real) = Float64(v)

Headers.headerval(::MolecularModelHeader, ::Val{:cell}, v::TriclinicPBC) = v

Headers.headerval(::MolecularModelHeader, ::Val{:pressure},
		v::AbstractMatrix{<:Real}) = Basis3D(v)

Headers.headerval(::MolecularModelHeader, ::Val{:virial},
		v::AbstractMatrix{<:Real}) = Basis3D(v)

mutable struct MolecularModel <: Multicollection
	header::MolecularModelHeader
	n::Int
	D::Dict{Symbol,Any}

	MolecularModel(n::Integer = 0) =
			new(MolecularModelHeader(), n, Dict{Symbol,Any}())
end

const MolecularModelView = MulticollectionView{MolecularModel}

const ParticleCollection = Union{MolecularModel,MolecularModelView}

const Particle = MulticollectionItem{MolecularModel}

abstract type MolecularTrajectory <: FormattedStream{MolecularModel} end

function MolecularModel(model::ParticleCollection)
	model2 = MolecularModel(length(model))
	for (key, val) in pairs(model)
		get!(model2, key, val)
	end
	merge!(model2.header, model.header)
	model2
end

MolecularModel(p::Particle) = MolecularModel(view(parent(p), parentindices(p)))

Base.similar(model::ParticleCollection, n::Integer) = MolecularModel(n)

function Multicollections.collval(::MolecularModel, ::Val{:topology},
		n::Integer, v::AbstractGraph)
	@boundscheck length(v) == n || error("expected $(n)-index graph")
	merge!(Graph(n), v)
end

Multicollections.collval(::MolecularModel, ::Val{:topology}, n::Integer,
		V::AbstractVector) = pair!(Graph(n), V...)

Multicollections.collval(::MolecularModel, ::Val{:ids}, n::Integer, v) =
		(Vector{Int}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:ids}, n::Integer,
		::UndefInitializer) = Vector{Int}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:names}, n::Integer, v) =
		(Vector{String}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:names}, n::Integer,
		::UndefInitializer) = Vector{String}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:resids}, n::Integer, v) =
		(Vector{Int}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:resids}, n::Integer,
		::UndefInitializer) = Vector{Int}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:resnames}, n::Integer, v) =
		(Vector{String}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:resnames}, n::Integer,
		::UndefInitializer) = Vector{String}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:chainids}, n::Integer, v) =
		(Vector{String}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:chainids}, n::Integer,
		::UndefInitializer) = Vector{String}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:elements}, n::Integer, v) =
		(Vector{String}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:elements}, n::Integer,
		::UndefInitializer) = Vector{String}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:R}, n::Integer, v) =
		(Vector{Vector3D}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:R}, n::Integer,
		::UndefInitializer) = Vector{Vector3D}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer, v) =
		(Vector{Vector3D}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer,
		::UndefInitializer) = Vector{Vector3D}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer, v) =
		(Vector{Vector3D}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer,
		::UndefInitializer) = Vector{Vector3D}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:masses}, n::Integer, v) =
		(Vector{Float64}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:masses}, n::Integer,
		::UndefInitializer) = Vector{Float64}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:charges}, n::Integer, v) =
		(Vector{Float64}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:charges}, n::Integer,
		::UndefInitializer) = Vector{Float64}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:bfactors}, n::Integer, v) =
		(Vector{Float64}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:bfactors}, n::Integer,
		::UndefInitializer) = Vector{Float64}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:occupancies}, n::Integer, v) =
		(Vector{Float64}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:occupancies}, n::Integer,
		::UndefInitializer) = Vector{Float64}(undef, n)

Multicollections.collval(::MolecularModel, ::Val{:SS}, n::Integer, v) =
		(Vector{String}(undef, n) .= v)

Multicollections.collval(::MolecularModel, ::Val{:SS}, n::Integer,
		::UndefInitializer) = Vector{String}(undef, n)

Multicollections.itemtocollprop(::Particle, ::Val{:id}) = :ids
Multicollections.colltoitemprop(::Particle, ::Val{:ids}) = :id

Multicollections.itemtocollprop(::Particle, ::Val{:name}) = :names
Multicollections.colltoitemprop(::Particle, ::Val{:names}) = :name

Multicollections.itemtocollprop(::Particle, ::Val{:resid}) = :resids
Multicollections.colltoitemprop(::Particle, ::Val{:resids}) = :resid

Multicollections.itemtocollprop(::Particle, ::Val{:resname}) = :resnames
Multicollections.colltoitemprop(::Particle, ::Val{:resnames}) = :resname

Multicollections.itemtocollprop(::Particle, ::Val{:chainid}) = :chainids
Multicollections.colltoitemprop(::Particle, ::Val{:chainids}) = :chainid

Multicollections.itemtocollprop(::Particle, ::Val{:element}) = :elements
Multicollections.colltoitemprop(::Particle, ::Val{:elements}) = :element

Multicollections.itemtocollprop(::Particle, ::Val{:mass}) = :masses
Multicollections.colltoitemprop(::Particle, ::Val{:masses}) = :mass

Multicollections.itemtocollprop(::Particle, ::Val{:charge}) = :charges
Multicollections.colltoitemprop(::Particle, ::Val{:charges}) = :charge

Multicollections.itemtocollprop(::Particle, ::Val{:bfactor}) = :bfactors
Multicollections.colltoitemprop(::Particle, ::Val{:bfactors}) = :bfactor

Multicollections.itemtocollprop(::Particle, ::Val{:occupancy}) = :occupancies
Multicollections.colltoitemprop(::Particle, ::Val{:occupancies}) = :occupancy

Multicollections.itemtocollprop(::Particle, ::Val{:ss}) = :SS
Multicollections.colltoitemprop(::Particle, ::Val{:SS}) = :ss

Multicollections.collitemname(::MolecularModel) = "particle"

function mcrptree(model::ParticleCollection)
	n = length(model)
	chainids = get(model, :chainids, Repeat("", n))
	resids = get(model, :resids, Repeat(0, n))
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

#=
for chain in chains(model)
	chain.chainid[1]
	for res in residues(chain)
		res.resid[1]
		res.resname[1]
		for p in res
			...
		end
	end
end

for res in residues(model)
	...
end

for frag in fragments(model)
	...
end

hierarchy = mcrp(model)
for (i, chain) in enumerate(chains(model, hierarchy))
	for (j, res) in enumerate(residues(model, i, hierarchy))
		for p in res
			...
		end
	end
end
=#

function chainat(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	chainids = get(model, :chainids, Repeat("", length(model)))
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
	chainids = get(model, :chainids, Repeat("", n))
	resids = get(model, :resids, Repeat(0, n))
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
