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

Headers.headerval(::MolecularModelHeader, ::Val{:cell},
        v::AbstractMatrix{<:Real}) = pbccell(v)

Headers.headerval(::MolecularModelHeader, ::Val{:pressure},
        v::AbstractMatrix{<:Real}) = pbccell(v)

Headers.headerval(::MolecularModelHeader, ::Val{:virial},
        v::AbstractMatrix{<:Real}) = pbccell(v)

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
        (FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:R}, n::Integer,
        ::UndefInitializer) =
        FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3))

Multicollections.collval(::MolecularModel, ::Val{:Rk}, n::Integer, v) =
        (FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:Rk}, n::Integer,
        ::UndefInitializer) =
        FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3))

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer, v) =
        (FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:V}, n::Integer,
        ::UndefInitializer) =
        FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3))

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer, v) =
        (FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3)) .= v)

Multicollections.collval(::MolecularModel, ::Val{:F}, n::Integer,
        ::UndefInitializer) =
        FixedArray(VectorBasedMatrix(Vector{Float64}(undef, 3*n), 3))

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

Multicollections.itemtocollprop(::Particle, ::Val{:R}) = :R
Multicollections.colltoitemprop(::Particle, ::Val{:R}) = :R

Multicollections.itemtocollprop(::Particle, ::Val{:Rk}) = :Rk
Multicollections.colltoitemprop(::Particle, ::Val{:Rk}) = :Rk

Multicollections.itemtocollprop(::Particle, ::Val{:V}) = :V
Multicollections.colltoitemprop(::Particle, ::Val{:V}) = :V

Multicollections.itemtocollprop(::Particle, ::Val{:F}) = :F
Multicollections.colltoitemprop(::Particle, ::Val{:F}) = :F

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
	chainids = get(model, :chainids, ScalarArray("", n))
	resids = get(model, :resids, ScalarArray(0, n))
	hierarchy(model, chainids, resids)
end

function hierarchy(model::ParticleCollection,
		chainids::AbstractVector{<:AbstractString},
		resids::AbstractVector{<:Integer})
	chainindices = similar(model, Int)
	resindices = similar(model, Int)
	localresindices = similar(model, Int)
	chainparticleindices = UnitRange{Int}[]
	resparticleindices = UnitRange{Int}[]
	localresparticleranges = Vector{UnitRange{Int}}[]
    n = length(model)
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
		chainindices[1] = 1
		resindices[1] = 1
		localresindices[1] = 1
		i = 2
		while i <= n
			thischainid = chainids[i]
			thisresid = resids[i]
			if thischainid != lastchainid
				push!(chainparticleindices, chainstart:chainend)
				push!(resparticleindices, resstart:resend)
				if thislocalresindex > length(localresparticleranges)
					push!(localresparticleranges, [resstart:resend])
				else
					push!(localresparticleranges[thislocalresindex],
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
					push!(resparticleindices, resstart:resend)
					if thislocalresindex > length(localresparticleranges)
						push!(localresparticleranges, [resstart:resend])
					else
						push!(localresparticleranges[thislocalresindex],
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
			chainindices[i] = thischainindex
			resindices[i] = thisresindex
			localresindices[i] = thislocalresindex
			i += 1
		end
		push!(chainparticleindices, chainstart:chainend)
		push!(resparticleindices, resstart:resend)
		if thislocalresindex > length(localresparticleranges)
			push!(localresparticleranges, [resstart:resend])
		else
			push!(localresparticleranges[thislocalresindex],
					resstart:resend)
		end
	end
	localresparticleindices = [RangeVector(i) for i in localresparticleranges]
	(chains = MulticollectionSplit(model, chainparticleindices,
            chainindices), residues = MulticollectionSplit(model,
            resparticleindices, resindices), localresidues =
			MulticollectionSplit(model, localresparticleindices,
            localresindices))
end

eachchain(model::ParticleCollection) = hierarchy(model).chains

eachresidue(model::ParticleCollection) = hierarchy(model).residues

eachlocalresidue(model::ParticleCollection) = hierarchy(model).localresidues

eachfragment(model::ParticleCollection) =
		eachfragment(model, get(model, :topology, Graph(length(model))))

function eachfragment(model::ParticleCollection, G::AbstractGraph)
	fragindices = similar(G, Int)
	fragparticleindices = Vector{Int}[]
	visited = falses(length(G))
	i = 1
	thisfragindex = 1
	while true
		I = sort!(connected!(Int[], G, i, visited))
		push!(fragparticleindices, I)
		fragindices[I] .= thisfragindex
		thisfragindex += 1
		i = findnext(!, visited, i+1)
		if i == nothing
			break
		end
	end
	MulticollectionSplit(model, fragparticleindices, fragindices)
end

function chainindices(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	chainids = get(model, :chainids, ScalarArray("", length(model)))
	chainindices(chainids, i)
end

function chainindices(chainids::AbstractArray{<:AbstractString}, i::Integer)
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

function residueindices(model::ParticleCollection, i::Integer)
	@boundscheck checkbounds(model, i)
	n = length(model)
	chainids = get(model, :chainids, ScalarArray("", n))
	resids = get(model, :resids, ScalarArray(0, n))
	residueindices(chainids, resids, i)
end

function residueindices(chainids::AbstractArray{<:AbstractString},
		resids::AbstractArray{<:Integer}, i::Integer)
	n = length(resids)
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

fragmentindices(model::ParticleCollection, i::Integer) =
		fragmentindices(get(model, :topology, Graph(length(model))), i)

fragmentindices(G::AbstractGraph, i::Integer) = sort!(connected(G, i))
