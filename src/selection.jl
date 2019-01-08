# Particle selection routines

const SelectionCache = Dict{Symbol,Any}

function emptyframe!(cache::SelectionCache)
	for key in [:pbcpos, :posgrids, :SS]
		delete!(cache, key)
	end
end

function getpbcpos!(cache::SelectionCache, R::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing})
	get!(cache, :pbcpos) do
		Rw, Rkw = get!(cache, :pbcposbuffer) do
			similar(R), similar(R)
		end
		pbcpos!(Rw, Rkw, R, cell)
	end
end

function getposgrid!(cache::SelectionCache, Rw::AbstractVector{Vector3D},
		Rkw::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		d::Real)
	grids = get!(cache, :posgrids) do
		Dict{Int,PositionGrid}()
	end
	ceild = ceil(Int, d)
	get!(grids, ceild) do
		buffers = get!(cache, :posgridbuffers) do
			Dict{Int,PositionGrid}()
		end
		buffer = get!(buffers, ceild) do
			posgrid(cell)
		end
		posgrid!(buffer, Rw, Rkw, cell, ceild)
	end
end

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

selection_pos_f(f) = checkcallable(f)

selection_pos_f(R::AbstractVector{<:Real}) = selection_pos_f(Vector3D(R))

selection_pos_f(R::Vector3D) = selection_pos_f([R])

selection_pos_f(R::AbstractVector{Vector3D}) = (model, cache) -> R

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

selection_vector3d_predicate(f) = checkcallable(f)

selection_vector3d_predicate(V::AbstractVector{<:Real}) =
		selection_vector3d_predicate(Vector3D(V))

selection_vector3d_predicate(V::Vector3D) = (W -> isapprox(W, V))

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
	I = s.f(eachindex(split.Igroup2items))
	selectindices!(work, I, split.Igroup2items)
	for i in subset
		results[i] = work[i]
	end
end

function selectindices!(results::BitVector, I::AbstractVector{<:Integer},
		Igroup2ps::AbstractVector{<:AbstractVector{<:Integer}})
	for i in I
		for j in Igroup2ps[i]
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
	selectpartial!(work, results, split.Igroup2items, split.Iitem2group, subset)
	for i in subset
		results[i] = work[i]
	end
end

function selectpartial!(results::BitVector, selected::BitVector,
		Igroup2ps::AbstractVector{<:AbstractVector{<:Integer}},
		Ip2group::AbstractVector{<:Integer}, subset::AbstractVector{<:Integer})
	for i in subset
		if selected[i]
			if ! results[i]
				for j in Igroup2ps[Ip2group[i]]
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
	selectfull!(work, results, split.Igroup2items, split.Iitem2group, subset)
	for i in subset
		results[i] = work[i]
	end
end

function selectfull!(results::BitVector, selected::BitVector,
		Igroup2ps::AbstractVector{<:AbstractVector{<:Integer}},
		Ip2group::AbstractVector{<:Integer},
		subset::AbstractVector{<:Integer})
	for i in subset
		if ! selected[i]
			if results[i]
				for j in Igroup2ps[Ip2group[i]]
					results[j] = false
				end
			end
		end
	end
end

struct WithinPositionSelector{T} <: Selector
	d::Float64
	of::T
	by::SelectBy

	function WithinPositionSelector{T}(d::Real, of::T, by::SelectBy) where {T}
		d > 0.0 || error("expected strictly positive distance")
		new(d, of, by)
	end
end

WithinPositionSelector(d::Real, of::T, by::SelectBy) where {T} =
		WithinPositionSelector{T}(d, of, by)

WithinSelector(d::Real, of, by::SelectBy) =
		WithinPositionSelector(d, selection_pos_f(of), by)

Expand(s::WithinPositionSelector; by::SelectBy) =
		WithinPositionSelector(s.d, s.of, by)

function select!(s::WithinPositionSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	work = trues(length(model))
	for i in subset
		work[i] = false
	end
	split = s.by(model, cache)
	cell = get(model.header, :cell, nothing)
	Rw, Rkw = getpbcpos!(cache, model.R, cell)
	pg = getposgrid!(cache, Rw, Rkw, cell, s.d)
	Rref = s.of(model, cache)
	selectwithinpos!(work, split.Igroup2items, split.Iitem2group, Rw, Rkw, cell,
			Rref, s.d, pg)
	for i in subset
		results[i] = work[i]
	end
end

function selectwithinpos!(results::BitVector,
		Igroup2ps::AbstractVector{<:AbstractVector{<:Integer}},
		Ip2group::AbstractVector{<:Integer}, Rw::AbstractVector{Vector3D},
		Rkw::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		Rref::AbstractVector{Vector3D}, d::Real, pg::PositionGrid)
	d2 = d^2
	I = Int[]
	for Rrefi in Rref
		Rrefiw, Rrefikw = pbcpos(Rrefi, cell)
		for i in findnear!(I, pg, Rrefi)
			if ! results[i]
				if sqmindist(Rw[i], Rkw[i], Rrefiw, Rrefikw, cell) <= d2
					for j in Igroup2ps[Ip2group[i]]
						results[j] = true
					end
				end
			end
		end
	end
end

struct WithinSelectionSelector{T<:Selector} <: Selector
	d::Float64
	of::T
	by::SelectBy

	function WithinSelectionSelector{T}(d::Real, of::T, by::SelectBy) where
			{T<:Selector}
		d > 0.0 || error("expected strictly positive distance")
		new(d, of, by)
	end
end

WithinSelectionSelector(d::Real, of::T, by::SelectBy) where {T<:Selector} =
		WithinSelectionSelector{T}(d, of, by)

WithinSelector(d::Real, of::Selector, by::SelectBy) =
		WithinSelectionSelector(d, of, by)

Within(d::Real; of) = WithinSelector(d, of, SelectByParticle())

Expand(s::WithinSelectionSelector; by::SelectBy) =
		WithinSelectionSelector(s.d, s.of, by)

function select!(s::WithinSelectionSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache)
	work = trues(length(model))
	for i in subset
		work[i] = false
	end
	split = s.by(model, cache)
	cell = get(model.header, :cell, nothing)
	Rw, Rkw = getpbcpos!(cache, model.R, cell)
	pg = getposgrid!(cache, Rw, Rkw, cell, s.d)
	Iref = findall(map(s.of, model, cache))
	selectwithinsel!(work, split.Igroup2items, split.Iitem2group, Rw, Rkw, cell,
			Iref, s.d, pg)
	for i in subset
		results[i] = work[i]
	end
end

function selectwithinsel!(results::BitVector,
		Igroup2ps::AbstractVector{<:AbstractVector{<:Integer}},
		Ip2group::AbstractVector{<:Integer}, Rw::AbstractVector{Vector3D},
		Rkw::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		Iref::AbstractVector{<:Integer}, d::Real, pg::PositionGrid)
	d2 = d^2
	J = Int[]
	for i in Iref
		for j in findnear!(J, pg, i)
			if ! results[j]
				if sqmindist(Rw[j], Rkw[j], Rw[i], Rkw[i], cell) <= d2
					for k in Igroup2ps[Ip2group[j]]
						results[k] = true
					end
				end
			end
		end
	end
end

struct PredicateSelector{T1,T2} <: Selector
	getter::T1
	predicate::T2

	PredicateSelector{T1,T2}(getter::T1, predicate::T2) where {T1,T2} =
			new(getter, predicate)
end

PredicateSelector(getter::T1, predicate::T2) where {T1,T2} =
		PredicateSelector{T1,T2}(getter, predicate)

select!(s::PredicateSelector, results::BitVector,
		model::ParticleCollection, subset::AbstractVector{<:Integer},
		cache::SelectionCache) = selectpredicated!(s.predicate, results,
		s.getter(model, cache), subset)

function selectpredicated!(f, results::BitVector, V::AbstractVector,
		subset::AbstractVector{<:Integer})
	for i in subset
		results[i] = f(V[i])
	end
end

Id(X...) = PredicateSelector((model, cache) -> model.ids,
		selection_int_predicate(X...))

Name(X...) = PredicateSelector((model, cache) -> model.names,
		selection_string_predicate(X...))

ResId(X...) = PredicateSelector((model, cache) -> model.resids,
		selection_int_predicate(X...))

ResName(X...) = PredicateSelector((model, cache) -> model.resnames,
		selection_string_predicate(X...))

ChainId(X...) = PredicateSelector((model, cache) ->
		get(model, :chainids, Repeat("", length(model))),
		selection_string_predicate(X...))

function Element(X...)
	f = function (model, cache)
		get(model, :elements) do
			get!(cache, :elements) do
				guesselements!(similar(model, String), model)
			end
		end
	end
	PredicateSelector(f, selection_string_predicate(X...))
end

R(x) = PredicateSelector((model, cache) -> model.R,
		selection_vector3d_predicate(x))

V(x) = PredicateSelector((model, cache) -> model.V,
		selection_vector3d_predicate(x))

F(x) = PredicateSelector((model, cache) -> model.F,
		selection_vector3d_predicate(x))

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
	PredicateSelector(f, selection_real_predicate(x))
end

Charge(x) = PredicateSelector((model, cache) -> model.charges,
		selection_real_predicate(x))

BFactor(x) = PredicateSelector((model, cache) -> model.bfactors,
		selection_real_predicate(x))

Occupancy(x) = PredicateSelector((model, cache) -> model.occupancies,
		selection_real_predicate(x))

function SS(X...)
	f = function (model, cache)
		get(model, :SS) do
			get!(cache, :SS) do
				guessss!(fill("", length(model)), model)
			end
		end
	end
	PredicateSelector(f, selection_string_predicate(X...))
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
