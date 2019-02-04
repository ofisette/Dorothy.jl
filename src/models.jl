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

Headers.headerval(::MolecularModelHeader, ::Val{:cell}, v::TriclinicPBC) =
		TriclinicCell(v)

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
