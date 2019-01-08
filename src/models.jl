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
	merge!(FixedGraph(Graph(n)), v)
end

function Multicollections.collval(::MolecularModel, ::Val{:topology},
		n::Integer, V::AbstractVector)
	G = Graph(n)
	pair!(G, [(i, j) for (i, j) in V]...)
end

Multicollections.collval(::MolecularModel, ::Val{:ids}, n::Integer, v) =
		(FixedArray(Vector{Int}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:ids}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Int}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:names}, n::Integer, v) =
		(FixedArray(Vector{String}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:names}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{String}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:resids}, n::Integer, v) =
		(FixedArray(Vector{Int}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:resids}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Int}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:resnames}, n::Integer, v) =
		(FixedArray(Vector{String}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:resnames}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{String}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:chainids}, n::Integer, v) =
		(FixedArray(Vector{String}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:chainids}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{String}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:elements}, n::Integer, v) =
		(FixedArray(Vector{String}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:elements}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{String}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:R}, n::Integer, v) =
		(FixedArray(Vector{Vector3D}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:R}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Vector3D}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer, v) =
		(FixedArray(Vector{Vector3D}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Vector3D}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer, v) =
		(FixedArray(Vector{Vector3D}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Vector3D}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:masses}, n::Integer, v) =
		(FixedArray(Vector{Float64}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:masses}, n::Integer,
		::UndefInitializer) =
		FixedArray(Vector{Float64}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:charges}, n::Integer, v) =
		(FixedArray(Vector{Float64}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:charges}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Float64}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:bfactors}, n::Integer, v) =
		(FixedArray(Vector{Float64}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:bfactors}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Float64}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:occupancies}, n::Integer, v) =
		(FixedArray(Vector{Float64}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:occupancies}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{Float64}(undef, n))

Multicollections.collval(::MolecularModel, ::Val{:SS}, n::Integer, v) =
		(FixedArray(Vector{String}(undef, n)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:SS}, n::Integer,
		::UndefInitializer) = FixedArray(Vector{String}(undef, n))

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

function hierarchy(model::ParticleCollection)
	n = length(model)
	chainids = get(model, :chainids, Repeat("", n))
	resids = get(model, :resids, Repeat(0, n))
	(Ichain2ps, Ip2chain), (Ires2ps, Ip2res), (Ilocalres2ps, Ip2localres) =
			hierarchy(chainids, resids)
	(chains = MulticollectionSplit(model, Ichain2ps, Ip2chain),
			residues = MulticollectionSplit(model, Ires2ps, Ip2res),
			localresidues = MulticollectionSplit(model, Ilocalres2ps,
			Ip2localres))
end

function hierarchy(chainids::AbstractVector{<:AbstractString},
		resids::AbstractVector{<:Integer})
	n = length(chainids)
	@boundscheck length(resids) == n ||
			error("size mismatch between chainid and resid arrays")
	Ip2chain = Vector{Int}(undef, n)
	Ip2res = Vector{Int}(undef, n)
	Ip2localres = Vector{Int}(undef, n)
	Ichain2ps = UnitRange{Int}[]
	Ires2ps = UnitRange{Int}[]
	localresranges = Vector{UnitRange{Int}}[]
	if n > 0
		thischainindex = 1
		thisresindex = 1
		thislocalresindex = 1
		chainstart = 1
		chainend = 1
		resstart = 1
		resend = 1
		lastchainid = chainids[1]
		lastresid = resids[1]
		Ip2chain[1] = 1
		Ip2res[1] = 1
		Ip2localres[1] = 1
		i = 2
		while i <= n
			thischainid = chainids[i]
			thisresid = resids[i]
			if thischainid != lastchainid
				push!(Ichain2ps, chainstart:chainend)
				push!(Ires2ps, resstart:resend)
				if thislocalresindex > length(localresranges)
					push!(localresranges, [resstart:resend])
				else
					push!(localresranges[thislocalresindex],
							resstart:resend)
				end
				thischainindex += 1
				thisresindex += 1
				thislocalresindex = 1
				chainstart = i
				chainend = i
				resstart = i
				resend = i
				lastchainid = thischainid
				lastresid = thisresid
			else
				chainend += 1
				if thisresid != lastresid
					push!(Ires2ps, resstart:resend)
					if thislocalresindex > length(localresranges)
						push!(localresranges, [resstart:resend])
					else
						push!(localresranges[thislocalresindex],
								resstart:resend)
					end
					thisresindex += 1
					thislocalresindex += 1
					resstart = i
					resend = i
					lastresid = thisresid
				else
					resend += 1
				end
			end
			Ip2chain[i] = thischainindex
			Ip2res[i] = thisresindex
			Ip2localres[i] = thislocalresindex
			i += 1
		end
		push!(Ichain2ps, chainstart:chainend)
		push!(Ires2ps, resstart:resend)
		if thislocalresindex > length(localresranges)
			push!(localresranges, [resstart:resend])
		else
			push!(localresranges[thislocalresindex], resstart:resend)
		end
	end
	Ilocalres2ps = [RangeVector(i) for i in localresranges]
	(Ichain2ps, Ip2chain), (Ires2ps, Ip2res), (Ilocalres2ps, Ip2localres)
end

eachchain(model::ParticleCollection) = hierarchy(model).chains

eachresidue(model::ParticleCollection) = hierarchy(model).residues

eachlocalresidue(model::ParticleCollection) = hierarchy(model).localresidues

function eachfragment(model::ParticleCollection)
	G = get(model, :topology, Graph(length(model)))
	(Ifrag2ps, Ip2frag) = eachfragment(G)
	MulticollectionSplit(model, Ifrag2ps, Ip2frag)
end

function eachfragment(G::AbstractGraph)
	Ip2frag = similar(G, Int)
	Ifrag2ps = Vector{Int}[]
	visited = falses(length(G))
	i = 1
	thisfragindex = 1
	while true
		I = sort!(connected!(Int[], G, i, visited))
		push!(Ifrag2ps, I)
		Ip2frag[I] .= thisfragindex
		thisfragindex += 1
		i = findnext(!, visited, i+1)
		if i == nothing
			break
		end
	end
	Ifrag2ps, Ip2frag
end

function chainat(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	chainids = get(model, :chainids, Repeat("", length(model)))
	view(model, chainat(chainids, i))
end

function chainat(chainids::AbstractArray{<:AbstractString}, i::Integer)
	@boundscheck checkbounds(chainids, i)
	n = length(chainids)
	thischainid = chainids[i]
	firsti = i
	firsti -= 1
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
