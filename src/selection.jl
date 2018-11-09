# Particle selection routines

const SelectionCache = Dict{Symbol,Any}

abstract type SelectBy end

struct SelectByChain <: SelectBy end

const Chain = SelectByChain()

(::SelectByChain)(model::ParticleCollection, cache::SelectionCache) =
		get!(cache, :hierarchy, hierarchy(model)).chains

struct SelectByResidue <: SelectBy end

const Residue = SelectByResidue()

(::SelectByResidue)(model::ParticleCollection, cache::SelectionCache) =
		get!(cache, :hierarchy, hierarchy(model)).residues

struct SelectByLocalResidue <: SelectBy end

const LocalResidue = SelectByLocalResidue()

(::SelectByLocalResidue)(model::ParticleCollection, cache::SelectionCache) =
		get!(cache, :hierarchy, hierarchy(model)).localresidues

struct SelectByFragment <: SelectBy end

const Fragment = SelectByFragment()

(::SelectByFragment)(model::ParticleCollection, cache::SelectionCache) =
		get!(cache, :fragments, eachfragment(model))

struct SelectByParticle <: SelectBy end

(::SelectByParticle)(model::ParticleCollection, cache::SelectionCache) =
		eachitem(model)

abstract type Selector end

Base.map(s::Selector, model::ParticleCollection, cache = SelectionCache()) =
		map!(s, BitVector(undef, length(model)), model, cache)

function Base.map!(s::Selector, results::BitVector, model::ParticleCollection,
		cache = SelectionCache())
	if ! isempty(model)
		select!(s, results, model, eachindex(model), cache)
	end
	results
end

Base.view(model::ParticleCollection, s::Selector, cache = SelectionCache()) =
		MulticollectionView(model, findall(map(s, model, cache)))

function pick(model::ParticleCollection, s::Selector, cache = SelectionCache())
	I = findall(map(s, model, cache))
	length(I) == 1 || error("expected exactly one particle")
	model[I[1]]
end

pick(model::ParticleCollection, chainid::AbstractString, resid::Real,
		name::AbstractString, cache = Selectioncache()) =
		pick(model, ChainId(chainid) & ResId(resid) & Name(name), cache)

pick(model::ParticleCollection, chainid::AbstractString, name::AbstractString,
		cache = SelectionCache()) =
		pick(model, ChainId(chainid) & Name(name), cache)

pick(model::ParticleCollection, resid::Real, name::AbstractString,
		cache = SelectionCache()) =
		pick(model, ResId(resid) & Name(name), cache)

pick(model::ParticleCollection, name::AbstractString,
		cache = SelectionCache()) = pick(model, Name(name), cache)

iscallable(f) = !isempty(methods(f))

function checkcallable(f)
	if ! iscallable(f)
		error("expected callable")
	end
	f
end

selection_index_f(f) = checkcallable(f)

selection_index_f(i::Integer) = (I -> [I[i]])

selection_index_f(I::AbstractVector{<:Integer}) = (J -> J[I])

selection_index_f(i0::Integer, I::Integer...) = selection_index_f([i0, I...])

selection_position_getter(R::AbstractMatrix{<:Real}) = (model, cache) -> R

selection_position_getter(R::AbstractVector{<:Real}) =
		(model, cache) -> reshape(R, 3,1)

selection_position_getter(s::Selector) =
		(model, cache) -> view(model, s, cache).R

selection_int_predicate(f) = checkcallable(f)

selection_int_predicate(i::Integer) = (==(i))

function selection_int_predicate(I::AbstractVector{<:Integer})
	J = sort(I)
	i -> !isempty(searchsorted(J, i))
end

selection_int_predicate(i0::Integer, I::Integer...) =
		selection_int_predicate([i0, I...])

selection_string_predicate(f) = checkcallable(f)

selection_string_predicate(s::AbstractString) = namematcher(s)

selection_string_predicate(s0::AbstractString, S::AbstractString...) =
		namematcher(s0, S...)

selection_real_predicate(f) = checkcallable(f)

selection_real_predicate(r::Real) = (s -> isapprox(s, r))

selection_vector_real_predicate(f) = checkcallable(f)

selection_vector_real_predicate(R::AbstractVector{<:Real}) =
		(S -> isapprox(R, S))

struct NotSelector{T<:Selector} <: Selector
	s::T

	NotSelector{T}(s::T) where {T<:Selector} = new(s)
end

function select!(s::NotSelector, results::BitVector, model::ParticleCollection,
		subset::AbstractVector{<:Integer}, cache::SelectionCache)
	select!(s.s, results, model, subset, cache)
	for i in subset
		results[i] = !results[i]
	end
end

Base.:(!)(s::T) where {T<:Selector} = NotSelector{T}(s)

struct AndSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::T1
	s2::T2

	AndSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(s::AndSelector, results::BitVector, model::ParticleCollection,
		subset::AbstractVector{<:Integer}, cache::SelectionCache)
	select!(s.s1, results, model, subset, cache)
	subset2 = [i for i in subset if results[i]]
	results2 = similar(results)
	select!(s.s2, results2, model, subset2, cache)
	for i in subset2
		if ! results2[i]
			results[i] = false
		end
	end
end

Base.:(&)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
		AndSelector{T1,T2}(s1, s2)

struct OrSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::T1
	s2::T2

	OrSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(s::OrSelector, results::BitVector, model::ParticleCollection,
		subset::AbstractVector{<:Integer}, cache::SelectionCache)
	select!(s.s1, results, model, subset, cache)
	subset2 = [i for i in subset if ! results[i]]
	results2 = similar(results)
	select!(s.s2, results2, model, subset2, cache)
	for i in subset2
		if results2[i]
			results[i] = true
		end
	end
end

Base.:(|)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
		OrSelector{T1,T2}(s1, s2)

struct XorSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::Selector
	s2::Selector

	XorSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(s::XorSelector, results::BitVector, model::ParticleCollection,
		subset::AbstractVector{<:Integer}, cache::SelectionCache)
	select!(s.s1, results, model, subset, cache)
	results2 = similar(results)
	select!(s.s2, results2, model, subset, cache)
	for i in subset
		results[i] = results[i] ⊻ results2[i]
	end
end

Base.:(⊻)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
		XorSelector{T1,T2}(s1, s2)

struct CachedSelector{T<:Selector} <: Selector
	s::T
	id::Symbol

	CachedSelector{T}(s::T, id::Symbol) where {T<:Selector} = new(s, id)
end

CachedSelector(s::T, id::Symbol) where {T<:Selector} = CachedSelector{T}(s, id)

function select!(s::CachedSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	cached::BitVector = get!(cache, s.id) do
		work = BitVector(undef, length(model))
		select!(s.s, work, model, eachindex(model), cache)
		work
	end
	for i in subset
		results[i] = cached[i]
	end
end

struct MapSelector{T<:AbstractVector{<:Bool}} <: Selector
	map::T

	MapSelector{T}(map::T) where {T<:AbstractVector{<:Bool}} = new(map)
end

MapSelector(map::T) where {T<:AbstractVector{<:Bool}} = MapSelector{T}(map)

Map(map::AbstractVector{<:Bool}) = MapSelector(map)

function select!(s::MapSelector, results::BitVector, model::ParticleCollection,
		subset::AbstractVector{<:Integer}, cache::SelectionCache)
	for i in subset
		results[i] = s.map[i]
	end
end

struct IndexSelector{T1,T2<:SelectBy} <:Selector
	f::T1
	by::T2

	IndexSelector{T1,T2}(f::T1, by::T2) where {T1,T2<:SelectBy} = new(f, by)
end

IndexSelector(f::T1, by::T2) where {T1,T2<:SelectBy} =
		IndexSelector{T1,T2}(f, by)

Index(X...; by::SelectBy = SelectByParticle()) =
		IndexSelector(selection_index_f(X...), by)

function select!(s::IndexSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	work = falses(length(model))
	split = s.by(model, cache)
	I = s.f(eachindex(split.itemindices))
	selectindices!(work, I, split.itemindices)
	for i in subset
		results[i] = work[i]
	end
end

function selectindices!(results::BitVector, I::AbstractVector{<:Integer},
		itemindices::AbstractVector{<:AbstractVector{<:Integer}})
	for i in I
		for j in itemindices[i]
			results[j] = true
		end
	end
end

struct ExpandSelector{T1<:Selector,T2<:SelectBy} <: Selector
	s::T1
	by::T2

	ExpandSelector{T1,T2}(s::T1, by::T2) where {T1<:Selector,T2<:SelectBy} =
			new(s, by)
end

ExpandSelector(s::T1, by::T2) where {T1<:Selector,T2<:SelectBy} =
		ExpandSelector{T1,T2}(s, by)

Expand(s::Selector; by::SelectBy) = ExpandSelector(s, by)

function select!(s::ExpandSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache,)
	select!(s.s, results, model, subset, cache)
	work = falses(length(model))
	split = s.by(model, cache)
	selectpartial!(work, results, split.itemindices, split.groupindices, subset)
	for i in subset
		results[i] = work[i]
	end
end

function selectpartial!(results::BitVector, selected::BitVector,
		itemindices::AbstractVector{<:AbstractVector{<:Integer}},
		groupindices::AbstractVector{<:Integer}, subset::AbstractVector{<:Integer})
	for i in subset
		if selected[i]
			if ! results[i]
				for j in itemindices[groupindices[i]]
					results[j] = true
				end
			end
		end
	end
end

struct RestrictSelector{T1<:Selector,T2<:SelectBy} <: Selector
	s::T1
	by::T2

	RestrictSelector{T1,T2}(s::T1, by::T2) where {T1<:Selector,T2<:SelectBy} =
			new(s, by)
end

RestrictSelector(s::T1, by::T2) where {T1<:Selector,T2<:SelectBy} =
		RestrictSelector{T1,T2}(s, by)

Restrict(s::Selector; by::SelectBy) = RestrictSelector(s, by)

function select!(s::RestrictSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	select!(s.s, results, model, subset, cache)
	work = trues(length(model))
	split = s.by(model, cache)
	selectfull!(work, results, split.itemindices, split.groupindices, subset)
	for i in subset
		results[i] = work[i]
	end
end

function selectfull!(results::BitVector, selected::BitVector,
		itemindices::AbstractVector{<:AbstractVector{<:Integer}},
		groupindices::AbstractVector{<:Integer},
		subset::AbstractVector{<:Integer})
	for i in subset
		if ! selected[i]
			if results[i]
				for j in itemindices[groupindices[i]]
					results[j] = false
				end
			end
		end
	end
end

struct WithinSelector{T} <: Selector
	d::Float64
	getter!::T
	by::SelectBy
	pbc::Bool

	function WithinSelector{T}(d::Real, getter!::T, by::SelectBy,
			pbc::Bool) where {T}
		d > 0.0 || error("expected strictly positive distance")
		new(d, getter!, by, pbc)
	end
end

WithinSelector(d::Real, getter!::T, by::SelectBy, pbc::Bool) where {T} =
		WithinSelector{T}(d, getter!, by, pbc)

Within(d::Real, x; by::SelectBy = SelectByParticle(), pbc::Bool = true) =
		WithinSelector(d, selection_position_getter(x), by, pbc)

function select!(s::WithinSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	work = trues(length(model))
	for i in subset
		work[i] = false
	end
	S = s.getter!(model, cache)
	split = s.by(model, cache)
	cell = get(model.header, :cell, nothing)
	if s.pbc && cell != nothing
		Rk = get(model, :Rk) do
			get!(cache, :Rk) do
				wrappos!(inv(cell) * model.R)
			end
		end
		selectwithin!(work, split.itemindices, split.groupindices, Rk, S, s.d,
				cell, cache)
	else
		selectwithin!(work, split.itemindices, split.groupindices, model.R, S,
				s.d, cache)
	end
	for i in subset
		results[i] = work[i]
	end
end

function selectwithin!(results::BitVector,
		itemindices::AbstractVector{<:AbstractVector{<:Integer}},
		groupindices::AbstractVector{<:Integer}, R::AbstractMatrix{<:Real},
		S::AbstractMatrix{<:Real}, d::Real, cache::SelectionCache)
	spatial::Dict{Symbol,Any} = get!(cache, :spatial, Dict{Symbol,Any}())
	grids::Dict{Int,OrthoGrid3{Int}} =
			get!(spatial, :grids, Dict{Int,OrthoGrid3{Int}}())
	dgrid = ceil(Int, d)
	g3::OrthoGrid3{Int} = get!(grids, dgrid) do
		append!(OrthoGrid3{Int}(extent(R), dgrid), R, 1:ncols(R))
	end
	d2 = d^2
	J = Int[]
	for i = 1:ncols(S)
		B = @view S[:,i]
		for j in findnear!(J, g3, B)
			if ! results[j]
				A = @view R[:,j]
				if sqdist(A, B) <= d2
					for k in itemindices[groupindices[j]]
						results[k] = true
					end
				end
			end
		end
	end
end

function selectwithin!(results::BitVector,
		itemindices::AbstractVector{<:AbstractVector{<:Integer}},
		groupindices::AbstractVector{<:Integer}, Rk::AbstractMatrix{<:Real},
		S::AbstractMatrix{<:Real}, d::Real, cell::AbstractMatrix{<:Real},
		cache::SelectionCache)
	spatial::Dict{Symbol,Any} = get!(cache, :spatial, Dict{Symbol,Any}())
	kgrids::Dict{Int,KspaceGrid3{Int}} =
			get!(spatial, :kgrids, Dict{Int,KspaceGrid3{Int}}())
	dgrid = ceil(Int, d)
	g3::KspaceGrid3{Int} = get!(kgrids, dgrid) do
		append!(KspaceGrid3{Int}(cell, dgrid), Rk, 1:ncols(Rk))
	end
	d2 = d^2
	Bk = similar(Rk, 3)
	Tk = similar(Rk, 3)
	T = similar(Rk, 3)
	J = Int[]
	for i = 1:ncols(S)
		B = @view S[:,i]
		mul!(Bk, cell.inv, B)
		wrappos!(Bk)
		for j in findnear!(J, g3, Bk)
			if ! results[j]
				Ak = @view Rk[:,j]
				mintrans!(Tk, Ak, Bk)
				mul!(T, cell, Tk)
				if sqnorm(T) <= d2
					for k in itemindices[groupindices[j]]
						results[k] = true
					end
				end
			end
		end
	end
end

struct VectorPredicateSelector{T1,T2} <: Selector
	getter!::T1
	predicate::T2

	VectorPredicateSelector{T1,T2}(getter!::T1, predicate::T2) where {T1,T2} =
			new(getter!, predicate)
end

VectorPredicateSelector(getter!::T1, predicate::T2) where{T1,T2} =
		VectorPredicateSelector{T1,T2}(getter!, predicate)

select!(s::VectorPredicateSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache) = selectpredicated!(s.predicate, results,
		s.getter!(model, cache), subset)

function selectpredicated!(f, results::BitVector, V::AbstractVector,
		subset::AbstractVector{<:Integer})
	for i in subset
		results[i] = f(V[i])
	end
end

struct MatrixPredicateSelector{T1,T2} <: Selector
	getter!::T1
	predicate::T2

	MatrixPredicateSelector{T1,T2}(getter!::T1, predicate::T2) where {T1,T2} =
			new(getter!, predicate)
end

MatrixPredicateSelector(getter!::T1, predicate::T2) where {T1,T2} =
		MatrixPredicateSelector{T1,T2}(getter!, predicate)

select!(s::MatrixPredicateSelector, results::BitVector,
		subset::AbstractVector{<:Integer}, model::ParticleCollection,
		cache::SelectionCache) = selectpredicated!(s.predicate, results,
		s.getter!(model, cache), subset)

function selectpredicated!(f, results::BitVector, M::AbstractMatrix,
		subset::AbstractVector{<:Integer})
	for i in subset
		results[i] = f(view(M, :, i))
	end
end

Id(X...) = VectorPredicateSelector((model, cache) -> model.ids,
		selection_int_predicate(X...))

Name(X...) = VectorPredicateSelector((model, cache) -> model.names,
		selection_string_predicate(X...))

ResId(X...) = VectorPredicateSelector((model, cache) -> model.resids,
		selection_int_predicate(X...))

ResName(X...) = VectorPredicateSelector((model, cache) -> model.resnames,
		selection_string_predicate(X...))

ChainId(X...) = VectorPredicateSelector((model, cache) ->
		get(model, :chainids, ScalarArray("", length(model))),
		selection_string_predicate(X...))

function Element(X...)
	f = function (model, cache)
		get(model, :elements) do
			get!(cache, :elements) do
				guesselements!(similar(model, String), model)
			end
		end
	end
	VectorPredicateSelector(f, selection_string_predicate(X...))
end

R(x) = MatrixPredicateSelector((model, cache) -> model.R,
		selection_vector_real_predicate(x))

V(x) = MatrixPredicateSelector((model, cache) -> model.V,
		selection_vector_real_predicate(x))

F(x) = MatrixPredicateSelector((model, cache) -> model.F,
		selection_vector_real_predicate(x))

function Mass(x)
	f = function (model, cache)
		get(model, :masses) do
			elements = get!(cache, :elements) do
				guesselements!(similar(model, String), model)
			end
			get!(cache, :masses) do
				guessmasses!(similar(model, Float64), model.names, elements)
			end
		end
	end
	VectorPredicateSelector(f, selection_real_predicate(x))
end

Charge(x) = VectorPredicateSelector((model, cache) -> model.charges,
		selection_real_predicate(x))

BFactor(x) = VectorPredicateSelector((model, cache) -> model.bfactors,
		selection_real_predicate(x))

Occupancy(x) = VectorPredicateSelector((model, cache) -> model.occupancies,
		selection_real_predicate(x))

function SS(X...)
	f = function (model, cache)
		get(model, :SS) do
			get!(cache, :SS) do
				guessss!(fill("", length(model)), model)
			end
		end
	end
	VectorPredicateSelector(f, selection_string_predicate(X...))
end

const Hydrogen = CachedSelector(Name(ishydrogen), :hydrogens)

const Heavy = !Hydrogen

const VSite = CachedSelector(Name(isvsite), :vsites)

const Water = CachedSelector(ResName(iswater), :water)

const Protein = CachedSelector(ResName(isprotein), :protein)

const AcidResidue = CachedSelector(ResName(isacidresidue), :acidresidues)

const BasicResidue = CachedSelector(ResName(isbasicresidue), :basicresidues)

const ChargedResidue =
		CachedSelector(ResName(ischargedresidue), :chargedresidues)

const PolarResidue = CachedSelector(ResName(ispolarresidue), :polarresidues)

const HydrophobicResidue =
		CachedSelector(ResName(ishydrophobicresidue), :hydrophobicresidues)

const MainChainName = CachedSelector(Name(ismainchainname), :mainchainnames)

const MainChain = MainChainName & Protein

const SideChain = !MainChainName & Protein

const BackboneName = CachedSelector(Name(isbackbonename), :backbonenames)

const Backbone = BackboneName & Protein

const CAlpha = Name("CA") & Protein

const Cα = CAlpha

const Nter = Index(1, by=LocalResidue) & Protein

const Cter = Index(last, by=LocalResidue) & Protein

const NuclAcid = CachedSelector(ResName(isnuclacid), :nuclacid)

const Lipid = CachedSelector(ResName(islipid), :lipid)

const Ion = CachedSelector(ResName(ision), :ion)

const Helix = CachedSelector(SS(ishelix), :helix)

const AlphaHelix = SS(isalphahelix)

const Helix310 = SS(ishelix310)

const PiHelix = SS(ispihelix)

const Turn = SS(isturn)

const Sheet = CachedSelector(SS(issheet), :sheet)

const Strand = SS(isstrand)

const Bridge = SS(isbridge)

const Loop = CachedSelector(SS(isloop), :loop)

const Coil = SS(iscoil)

const Bend = SS(isbend)

module Selectors

	using ..Dorothy

	export
			Chain, Residue, LocalResidue, Fragment,

			Map, Index, Expand, Restrict, Within,

			Id, Name, ResId, ResName, ChainId, Element, R, V, F, Mass, Charge,
			BFactor, Occupancy, Charge, SS,

			Hydrogen, Heavy, VSite, Water, Protein, AcidResidue, BasicResidue,
			ChargedResidue, PolarResidue, HydrophobicResidue, MainChain,
			SideChain, Backbone, Nter, Cter, NuclAcid, Lipid, Ion, Helix,
			AlphaHelix, Helix310, PiHelix, Turn, Sheet, Strand, Bridge, Loop,
			Coil, Bend

	const Chain = Dorothy.Chain
	const Residue = Dorothy.Residue
	const LocalResidue = Dorothy.LocalResidue
	const Fragment = Dorothy.Fragment
	const Map = Dorothy.Map
	const Index = Dorothy.Index
	const Expand = Dorothy.Expand
	const Restrict = Dorothy.Restrict
	const Within = Dorothy.Within
	const Id = Dorothy.Id
	const Name = Dorothy.Name
	const ResId = Dorothy.ResId
	const ResName = Dorothy.ResName
	const ChainId = Dorothy.ChainId
	const Element = Dorothy.Element
	const R = Dorothy.R
	const V = Dorothy.V
	const F = Dorothy.F
	const Mass = Dorothy.Mass
	const Charge = Dorothy.Charge
	const BFactor = Dorothy.BFactor
	const Occupancy = Dorothy.Occupancy
	const Charge = Dorothy.Charge
	const SS = Dorothy.SS
	const Hydrogen = Dorothy.Hydrogen
	const Heavy = Dorothy.Heavy
	const VSite = Dorothy.VSite
	const Water = Dorothy.Water
	const Protein = Dorothy.Protein
	const AcidResidue = Dorothy.AcidResidue
	const BasicResidue = Dorothy.BasicResidue
	const ChargedResidue = Dorothy.ChargedResidue
	const PolarResidue = Dorothy.PolarResidue
	const HydrophobicResidue = Dorothy.HydrophobicResidue
	const MainChain = Dorothy.MainChain
	const SideChain = Dorothy.SideChain
	const Backbone = Dorothy.Backbone
	const CAlpha = Dorothy.CAlpha
	const Cα = Dorothy.Cα
	const Nter = Dorothy.Nter
	const Cter = Dorothy.Cter
	const NuclAcid = Dorothy.NuclAcid
	const Lipid = Dorothy.Lipid
	const Ion = Dorothy.Ion
	const Helix = Dorothy.Helix
	const AlphaHelix = Dorothy.AlphaHelix
	const Helix310 = Dorothy.Helix310
	const PiHelix = Dorothy.PiHelix
	const Turn = Dorothy.Turn
	const Sheet = Dorothy.Sheet
	const Strand = Dorothy.Strand
	const Bridge = Dorothy.Bridge
	const Loop = Dorothy.Loop
	const Coil = Dorothy.Coil
	const Bend = Dorothy.Bend

end
