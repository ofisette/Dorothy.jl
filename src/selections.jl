# Particle selection routines

struct SelectionCache
	D::Dict{Symbol,Any}
end

SelectionCache() = SelectionCache(Dict{Symbol,Any}())

function Base.get!(f, cache::SelectionCache, groupkey::Symbol, itemkey::Symbol)
	group = get!(cache.D, groupkey) do
		Dict{Symbol,Any}()
	end
	get!(f, group, itemkey)
end

Base.delete!(cache::SelectionCache, groupkey::Symbol) =
		delete!(cache.D, groupkey)

function getchainids(model::ParticleCollection, cache::SelectionCache)
	get(model, :chainids) do
		Repeated("", length(model))
	end
end

function getelements(model::ParticleCollection, cache::SelectionCache)
	get!(cache, :system, :elements) do
		infermissingelements(model)
	end
end

function getmasses(model::ParticleCollection, cache::SelectionCache)
	get(model, :masses) do
		get!(cache, :system, :masses) do
			masses = Vector{Float64}(undef, length(model))
			elements = getelements(model, cache)
			infermasses!(masses, model.names, elements)
		end
	end
end

function getss(model::ParticleCollection, cache::SelectionCache)
	get(model, :SS) do
		get!(cache, :frame, :SS) do
			inferss(model)
		end
	end
end

function getpbccell(model::ParticleCollection, cache::SelectionCache)
	get!(cache, :frame, :cell) do
		pbccell(get(model.header, :cell, nothing))
	end
end

function getpbcpos(R::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, cache::SelectionCache)
	get!(cache, :frame, :pbcpos) do
		Rw, Kw = get!(cache, :buffers, :pbcpos) do
			similar(R), similar(R)
		end
		pbcpos!(Rw, Kw, R, cell)
	end
end

function getnbrtable(n::Integer, cell::Union{TriclinicPBC,Nothing},
		cache::SelectionCache)
	get!(cache, :frame, :nbrtable) do
		buffer = get!(cache, :buffers, :nbrtable) do
			NeighborTableBuffer()
		end
		NeighborTable(n, cell, buffer)
	end
end

abstract type Selector end

Base.map(s::Selector, model::ParticleCollection, cache = SelectionCache()) =
		map!(s, BitVector(undef, length(model)), model, cache)

function Base.map!(s::Selector, results::BitVector, model::ParticleCollection,
		cache = SelectionCache())
	if ! isempty(model)
		select!(results, eachindex(model), s, model, cache)
	end
	results
end

Base.view(model::ParticleCollection, s::Selector, cache = SelectionCache()) =
		MulticollectionView(model, findall(map(s, model, cache)))

Base.getindex(model::ParticleCollection, s::Selector) = view(model, s)

abstract type SelectionMode end

struct ModelMode <: SelectionMode end

Base.show(io::IO, ::ModelMode) = "model"

struct ChainMode <: SelectionMode end

Base.show(io::IO, ::ChainMode) = "chain"

function getmcrp(model::ParticleCollection, cache::SelectionCache)
	get!(cache, :system, :mcrp) do
		mcrp(model)
	end
end

geth2(model::ParticleCollection, ::ChainMode, cache::SelectionCache) =
		getmcrp(model, cache).flath1

struct ResidueMode <: SelectionMode end

Base.show(io::IO, ::ResidueMode) = "residue"

geth2(model::ParticleCollection, ::ResidueMode, cache::SelectionCache) =
		getmcrp(model, cache).flath2

struct ParticleMode <: SelectionMode end

Base.show(io::IO, ::ParticleMode) = "particle"

struct FragmentMode <: SelectionMode end

Base.show(io::IO, ::FragmentMode) = "fragment"

function gettopology(model::ParticleCollection, cache::SelectionCache)
	get(model, :topology) do
		get!(cache, :system, :topology) do
			infertopology(model)
		end
	end
end

function geth2(model::ParticleCollection, ::FragmentMode, cache::SelectionCache)
	get!(cache, :system, :mfp) do
		mfp(mfptree(gettopology(model, cache)), length(model))
	end
end

iscallable(f) = !isempty(methods(f))

function checkcallable(f)
	if ! iscallable(f)
		error("expected callable")
	end
	f
end

indexfunction(f) = checkcallable(f)

indexfunction(i::Integer) = (J -> i)

indexfunction(I::AbstractVector{<:Integer}) = (J -> I)

positionfunction(f) = checkcallable(f)

positionfunction(R::AbstractVector{<:Real}) = positionfunction(Vector3D(R))

positionfunction(R::Vector3D) = positionfunction([R])

positionfunction(R::AbstractVector{Vector3D}) = (model, cache) -> R

integerpredicate(f) = checkcallable(f)

integerpredicate(i::Integer) = (==(i))

function integerpredicate(I::AbstractVector{<:Integer})
	J = sort(I)
	i -> !isempty(searchsorted(J, i))
end

stringpredicate(f) = checkcallable(f)

stringpredicate(s::AbstractString) = namematcher(s)

stringpredicate(S::AbstractVector{<:AbstractString}) = namematcher(S)

realpredicate(f) = checkcallable(f)

realpredicate(r::Real) = (s -> isapprox(s, r))

vector3dpredicate(f) = checkcallable(f)

vector3dpredicate(V::AbstractVector{<:Real}) = vector3dpredicate(Vector3D(V))

vector3dpredicate(V::Vector3D) = (W -> isapprox(W, V))

struct NotSelector{T<:Selector} <: Selector
	s::T

	NotSelector{T}(s::T) where {T<:Selector} = new(s)
end

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::NotSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s, model, cache)
	for i in subset
		results[i] = !results[i]
	end
end

struct AndSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::T1
	s2::T2

	AndSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::AndSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s1, model, cache)
	subset2 = [i for i in subset if results[i]]
	results2 = similar(results)
	select!(results2, subset2, s.s2, model, cache)
	for i in subset2
		if ! results2[i]
			results[i] = false
		end
	end
end

struct OrSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::T1
	s2::T2

	OrSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::OrSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s1, model, cache)
	subset2 = [i for i in subset if ! results[i]]
	results2 = similar(results)
	select!(results2, subset2, s.s2, model, cache)
	for i in subset2
		if results2[i]
			results[i] = true
		end
	end
end

struct XorSelector{T1<:Selector,T2<:Selector} <: Selector
	s1::Selector
	s2::Selector

	XorSelector{T1,T2}(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			new(s1, s2)
end

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::XorSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s1, model, cache)
	results2 = similar(results)
	select!(results2, subset2, s.s2, model, cache)
	for i in subset
		results[i] = results[i] ⊻ results2[i]
	end
end

struct CachedSelector{T<:Selector} <: Selector
	s::T
	groupkey::Symbol
	itemkey::Symbol

	CachedSelector{T}(s::T, groupkey::Symbol, itemkey::Symbol) where
			{T<:Selector} = new(s, groupkey, itemkey)
end

CachedSelector(s::T, groupkey::Symbol, itemkey::Symbol) where {T<:Selector} =
		CachedSelector{T}(s, groupkey, itemkey)

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::CachedSelector, model::ParticleCollection, cache::SelectionCache)
	cached::BitVector = get!(cache, s.groupkey, s.itemkey) do
		work = BitVector(undef, length(model))
		select!(work, eachindex(model), s.s, model, cache)
		work
	end
	for i in subset
		results[i] = cached[i]
	end
end

struct MapSelector{T<:AbstractVector{Bool}} <: Selector
	map::T

	MapSelector{T}(map::T) where {T<:AbstractVector{Bool}} = new(map)
end

MapSelector(map::T) where {T<:AbstractVector{Bool}} = MapSelector{T}(map)

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::MapSelector, model::ParticleCollection, cache::SelectionCache)
	for i in subset
		results[i] = s.map[i]
	end
end

struct IndexSelector{T1,T2<:SelectionMode,T3<:SelectionMode} <:Selector
	f::T1
	by::T2
	ineach::T3

	function IndexSelector{T1,T2,T3}(f::T1, by::T2, ineach::T3) where
			{T1,T2<:SelectionMode,T3<:SelectionMode}
		if ! iscompatible(IndexSelector, by, ineach)
			error("cannot select indices by $(by) in each $(ineach)")
		end
		new(f, by, ineach)
	end
end

IndexSelector(f::T1, by::T2, ineach::T3) where
		{T1,T2<:SelectionMode,T3<:SelectionMode} =
		IndexSelector{T1,T2,T3}(f, by, ineach)

iscompatible(::Type{<:IndexSelector},
		by::SelectionMode, ineach::SelectionMode) = false

iscompatible(::Type{<:IndexSelector},
		by::ParticleMode, ineach::ModelMode) = true

iscompatible(::Type{<:IndexSelector},
		by::SelectionMode, ineach::ModelMode) = true

iscompatible(::Type{<:IndexSelector},
		by::ParticleMode, ineach::SelectionMode) = true

iscompatible(::Type{<:IndexSelector}, by::ResidueMode, ineach::ChainMode) = true

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::IndexSelector, model::ParticleCollection, cache::SelectionCache)
	work = falses(length(model))
	selectindices!(work, s.f, s.by, s.ineach, model, cache)
	for i in subset
		results[i] = work[i]
	end
end

selectindices!(results::BitVector, f, by::ParticleMode, ineach::ModelMode,
		model::ParticleCollection, cache::SelectionCache) =
		selectindices0r!(results, f(eachindex(model)))

function selectindices0r!(results::BitVector,
		I::Union{AbstractVector{<:Integer},Integer})
	for i in I
		results[i] = true
	end
end

selectindices!(results::BitVector, f, by::ParticleMode, ineach::SelectionMode,
		model::ParticleCollection, cache::SelectionCache) =
		selectindices01!(results, f, geth2(model, ineach, cache).tree)

function selectindices01!(results::BitVector, f,
		tree::AbstractVector{<:AbstractVector{<:Integer}})
	for branch in tree
		selectindices01!(results, branch, f(eachindex(branch)))
	end
end

function selectindices01!(results::BitVector, branch::AbstractVector{<:Integer},
		I::Union{AbstractVector{<:Integer},Integer})
	for i in I
		results[branch[i]] = true
	end
end

selectindices!(results::BitVector, f, by::SelectionMode, ineach::ModelMode,
		model::ParticleCollection, cache::SelectionCache) =
		selectindices1r!(results, f, geth2(model, by, cache).tree)

selectindices1r!(results::BitVector, f,
		tree::AbstractVector{<:AbstractVector{<:Integer}}) =
		selectindices1r!(results, tree, f(eachindex(tree)))

function selectindices1r!(results::BitVector,
		tree::AbstractVector{<:AbstractVector{<:Integer}},
		I::Union{AbstractVector{<:Integer},Integer})
	for i in I
		for j in tree[i]
			results[j] = true
		end
	end
end

selectindices!(results::BitVector, f, by::ResidueMode, ineach::ChainMode,
		model::ParticleCollection, cache::SelectionCache) =
		selectindices12!(results, f, getmcrp(model, cache).tree)

function selectindices12!(results::BitVector, f,
		tree::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}})
	for branch in tree
		selectindices12!(results, branch, f(eachindex(branch)))
	end
end

function selectindices12!(results::BitVector,
		branch::AbstractVector{<:AbstractVector{<:Integer}},
		I::Union{AbstractVector{<:Integer},Integer})
	for i in I
		for j in branch[i]
			results[j] = true
		end
	end
end

abstract type ByGroupSelector <: Selector end

struct ExpandSelector{T1<:Selector,T2<:SelectionMode} <: ByGroupSelector
	s::T1
	by::T2

	function ExpandSelector{T1,T2}(s::T1, by::T2) where
			{T1<:Selector,T2<:SelectionMode}
		if ! iscompatible(ByGroupSelector, by)
			error("cannot expand selections by $(by)")
		end
		new(s, by)
	end
end

ExpandSelector(s::T1, by::T2) where {T1<:Selector,T2<:SelectionMode} =
		ExpandSelector{T1,T2}(s, by)

iscompatible(::Type{<:ByGroupSelector}, ::SelectionMode) = false

iscompatible(::Type{<:ByGroupSelector}, ::ChainMode) = true

iscompatible(::Type{<:ByGroupSelector}, ::ResidueMode) = true

iscompatible(::Type{<:ByGroupSelector}, ::FragmentMode) = true

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::ExpandSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s, model, cache)
	work = falses(length(model))
	selectpartial!(work, subset, results, geth2(model, s.by, cache)...)
	for i in subset
		results[i] = work[i]
	end
end

function selectpartial!(results::BitVector, subset::AbstractVector{<:Integer},
		selected::BitVector, tree::AbstractVector{<:AbstractVector{<:Integer}},
		paths::AbstractVector{<:Integer})
	for i in subset
		if selected[i]
			if ! results[i]
				for j in tree[paths[i]]
					results[j] = true
				end
			end
		end
	end
end

struct RestrictSelector{T1<:Selector,T2<:SelectionMode} <: ByGroupSelector
	s::T1
	by::T2

	function RestrictSelector{T1,T2}(s::T1, by::T2) where
			{T1<:Selector,T2<:SelectionMode}
			if ! iscompatible(ByGroupSelector, by)
				error("cannot restrict selections by $(by)")
			end
		new(s, by)
	end
end

RestrictSelector(s::T1, by::T2) where {T1<:Selector,T2<:SelectionMode} =
		RestrictSelector{T1,T2}(s, by)

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::RestrictSelector, model::ParticleCollection, cache::SelectionCache)
	select!(results, subset, s.s, model, cache)
	work = trues(length(model))
	selectfull!(work, subset, results, geth2(model, s.by, cache)...)
	for i in subset
		results[i] = work[i]
	end
end

function selectfull!(results::BitVector, subset::AbstractVector{<:Integer},
		selected::BitVector, tree::AbstractVector{<:AbstractVector{<:Integer}},
		paths::AbstractVector{<:Integer})
	for i in subset
		if ! selected[i]
			if results[i]
				for j in tree[paths[i]]
					results[j] = false
				end
			end
		end
	end
end

abstract type WithinSelector <: Selector end

struct WithinPositionSelector{T1,T2<:SelectionMode} <: WithinSelector
	dmax::Float64
	of::T1
	by::T2

	function WithinPositionSelector{T1,T2}(dmax::Real, of::T1, by::T2) where
			{T1,T2<:SelectionMode}
		if ! iscompatible(WithinSelector, by)
			error("cannot select within a distance by $(by)")
		end
		dmax > 0.0 || error("expected strictly positive distance")
		new(dmax, of, by)
	end
end

iscompatible(::Type{<:WithinSelector}, by::SelectionMode) =
		(by == ParticleMode() || iscompatible(ByGroupSelector, by))

WithinPositionSelector(dmax::Real, of::T1, by::T2) where
		{T1,T2<:SelectionMode} = WithinPositionSelector{T1,T2}(dmax, of, by)

WithinSelector(dmax::Real, of, by::SelectionMode) =
		WithinPositionSelector(dmax, positionfunction(of), by)

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::WithinPositionSelector, model::ParticleCollection,
		cache::SelectionCache)
	work = trues(length(model))
	for i in subset
		work[i] = false
	end
	cell = getpbccell(model, cache)
	Rw, Kw = getpbcpos(model.R, cell, cache)
	nbrtable = getnbrtable(length(Rw), cell, cache)
	Rref = s.of(model, cache)
	selectwithinpos!(work, s.dmax, s.by, Rw, Kw, Rref, cell, nbrtable, model,
			cache)
	for i in subset
		results[i] = work[i]
	end
end

selectwithinpos!(results::BitVector, dmax::Real, by::ParticleMode,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable, model::ParticleCollection,
		cache::SelectionCache) =
		selectwithinpos0!(results, dmax, Rw, Kw, Rref, cell, nbrtable)

function selectwithinpos0!(results::BitVector, dmax::Real,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable)
	dmax2 = dmax^2
	J = Int[]
	for Rrefi in Rref
		Rrefiw, Krefiw = pbcpos(Rrefi, cell)
		for j in findnear!(J, Rw, Kw, Rrefiw, Krefiw, cell, dmax, nbrtable)
			if ! results[j]
				if sqmindist(Rw[j], Kw[j], Rrefiw, Krefiw, cell) <= dmax2
					results[j] = true
				end
			end
		end
	end
end

selectwithinpos!(results::BitVector, dmax::Real, by::SelectionMode,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable, model::ParticleCollection,
		cache::SelectionCache) = selectwithinpos1!(results, dmax, Rw, Kw, Rref,
		cell, nbrtable, geth2(model, by, cache)...)

function selectwithinpos1!(results::BitVector, dmax::Real,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable,
		tree::AbstractVector{<:AbstractVector{<:Integer}},
		paths::AbstractVector{<:Integer})
	dmax2 = dmax^2
	J = Int[]
	for Rrefi in Rref
		Rrefiw, Krefiw = pbcpos(Rrefi, cell)
		for j in findnear!(J, Rw, Kw, Rrefiw, Krefiw, cell, dmax, nbrtable)
			if ! results[j]
				if sqmindist(Rw[j], Kw[j], Rrefiw, Krefiw, cell) <= dmax2
					for k in tree[paths[j]]
						results[k] = true
					end
				end
			end
		end
	end
end

struct WithinSelectionSelector{T1<:Selector,T2<:SelectionMode} <: WithinSelector
	dmax::Float64
	of::T1
	by::T2

	function WithinSelectionSelector{T1,T2}(dmax::Real, of::T1, by::T2) where
			{T1<:Selector,T2<:SelectionMode}
		if ! iscompatible(WithinSelector, by)
			error("cannot select within a distance by $(by)")
		end
		dmax > 0.0 || error("expected strictly positive distance")
		new(dmax, of, by)
	end
end

WithinSelectionSelector(dmax::Real, of::T1, by::T2) where
		{T1<:Selector,T2<:SelectionMode} =
		WithinSelectionSelector{T1,T2}(dmax, of, by)

WithinSelector(dmax::Real, of::Selector, by::SelectionMode) =
		WithinSelectionSelector(dmax, of, by)

function select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::WithinSelectionSelector, model::ParticleCollection,
		cache::SelectionCache)
	work = trues(length(model))
	for i in subset
		work[i] = false
	end
	cell = getpbccell(model, cache)
	Rw, Kw = getpbcpos(model.R, cell, cache)
	nbrtable = getnbrtable(length(Rw), cell, cache)
	Iref = findall(map(s.of, model, cache))
	selectwithinsel!(work, s.dmax, s.by, Rw, Kw, Iref, cell, nbrtable, model,
			cache)
	for i in subset
		results[i] = work[i]
	end
end

selectwithinsel!(results::BitVector, dmax::Real, by::ParticleMode,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Iref::AbstractVector{<:Integer}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable, model::ParticleCollection,
		cache::SelectionCache) =
		selectwithinsel0!(results, dmax, Rw, Kw, Iref, cell, nbrtable)

function selectwithinsel0!(results::BitVector, dmax::Real,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Iref::AbstractVector{<:Integer}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable)
	dmax2 = dmax^2
	J = Int[]
	for i in Iref
		for j in findnear!(J, Rw, Kw, i, cell, dmax, nbrtable)
			if ! results[j]
				if sqmindist(Rw[j], Kw[j], Rw[i], Kw[i], cell) <= dmax2
					results[j] = true
				end
			end
		end
	end
end

selectwithinsel!(results::BitVector, dmax::Real, by::SelectionMode,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Iref::AbstractVector{<:Integer}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable, model::ParticleCollection,
		cache::SelectionCache) = selectwithinsel1!(results, dmax, Rw, Kw, Iref,
		cell, nbrtable, geth2(model, by, cache)...)

function selectwithinsel1!(results::BitVector, dmax::Real,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Iref::AbstractVector{<:Integer}, cell::Union{TriclinicPBC,Nothing},
		nbrtable::NeighborTable,
		tree::AbstractVector{<:AbstractVector{<:Integer}},
		paths::AbstractVector{<:Integer})
	dmax2 = dmax^2
	J = Int[]
	for i in Iref
		for j in findnear!(J, Rw, Kw, i, cell, dmax, nbrtable)
			if ! results[j]
				if sqmindist(Rw[j], Kw[j], Rw[i], Kw[i], cell) <= dmax2
					for k in tree[paths[j]]
						results[k] = true
					end
				end
			end
		end
	end
end

struct PropertySelector{T1,T2} <: Selector
	getter::T1
	predicate::T2

	PropertySelector{T1,T2}(getter::T1, predicate::T2) where {T1,T2} =
			new(getter, predicate)
end

PropertySelector(getter::T1, predicate::T2) where {T1,T2} =
		PropertySelector{T1,T2}(getter, predicate)

select!(results::BitVector, subset::AbstractVector{<:Integer},
		s::PropertySelector, model::ParticleCollection, cache::SelectionCache) =
		selecttrue!(results, subset, s.predicate, s.getter(model, cache))

function selecttrue!(results::BitVector, subset::AbstractVector{<:Integer},
		predicate, V::AbstractVector)
	for i in subset
		results[i] = predicate(V[i])
	end
end

module Selectors

	using ..Dorothy
	using ..Dorothy.Utils
	using ..Dorothy.Properties

	export
			Chain, Residue, Fragment,

			Map, Index, Expand, Restrict, Within,

			Id, Name, ResId, ResName, ChainId, Element, R, V, F, Mass, Charge,
			BFactor, Occupancy, Charge, SS,

			VSite, Hydrogen, Heavy, Water, Protein, AcidResidue, BasicResidue,
			ChargedResidue, PolarResidue, HydrophobicResidue, Backbone,
			MainChain, SideChain, Calpha, Cα, Nter, Cter, NuclAcid, Lipid, Ion,
			Helix, AlphaHelix, Helix310, PiHelix, Turn, Sheet, Strand, Bridge,
			Loop, Coil, Bend

	const Model = Dorothy.ModelMode()
	const Chain = Dorothy.ChainMode()
	const Residue = Dorothy.ResidueMode()
	const Particle = Dorothy.ParticleMode()
	const Fragment = Dorothy.FragmentMode()

	Base.:(!)(s::T) where {T<:Selector} = Dorothy.NotSelector{T}(s)

	Base.:(&)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			Dorothy.AndSelector{T1,T2}(s1, s2)

	Base.:(|)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			Dorothy.OrSelector{T1,T2}(s1, s2)

	Base.:(⊻)(s1::T1, s2::T2) where {T1<:Selector,T2<:Selector} =
			Dorothy.XorSelector{T1,T2}(s1, s2)

	Map(map::AbstractVector{Bool}) = Dorothy.MapSelector(map)

	function Map(view::MolecularModelView)
		model = parent(view)
		map = falses(length(model))
		for i in parentindices(view)
			map[i] = true
		end
		Dorothy.MapSelector(map)
	end

	Index(X; by::SelectionMode = Particle, ineach::SelectionMode = Model) =
			Dorothy.IndexSelector(Dorothy.indexfunction(X), by, ineach)

	Expand(s::Selector; by::SelectionMode) = Dorothy.ExpandSelector(s, by)

	Restrict(s::Selector; by::SelectionMode) = Dorothy.RestrictSelector(s, by)

	Within(dmax::Real; of) = Dorothy.WithinSelector(dmax, of, Particle)

	Expand(s::Dorothy.WithinPositionSelector; by::SelectionMode) =
			Dorothy.WithinPositionSelector(s.dmax, s.of, by)

	Expand(s::Dorothy.WithinSelectionSelector; by::SelectionMode) =
			Dorothy.WithinSelectionSelector(s.dmax, s.of, by)

	Id(X...) = Dorothy.PropertySelector((model, cache) -> model.ids,
			Dorothy.integerpredicate(X...))

	Name(X...) = Dorothy.PropertySelector((model, cache) -> model.names,
			Dorothy.stringpredicate(X...))

	ResId(X...) = Dorothy.PropertySelector((model, cache) -> model.resids,
			Dorothy.integerpredicate(X...))

	ResName(X...) = Dorothy.PropertySelector((model, cache) -> model.resnames,
			Dorothy.stringpredicate(X...))

	ChainId(X...) = Dorothy.PropertySelector(Dorothy.getchainids,
			Dorothy.stringpredicate(X...))

	Element(X...) = Dorothy.PropertySelector(Dorothy.getelements,
			Dorothy.stringpredicate(X...))

	R(x) = Dorothy.PropertySelector((model, cache) -> model.R,
			Dorothy.vector3dpredicate(x))

	V(x) = Dorothy.PropertySelector((model, cache) -> model.V,
			Dorothy.vector3dpredicate(x))

	F(x) = Dorothy.PropertySelector((model, cache) -> model.F,
			Dorothy.vector3dpredicate(x))

	Mass(x) = Dorothy.PropertySelector(Dorothy.getmasses,
			Dorothy.realpredicate(x))

	Charge(x) = Dorothy.PropertySelector((model, cache) -> model.charges,
			Dorothy.realpredicate(x))

	BFactor(x) = Dorothy.PropertySelector((model, cache) -> model.bfactors,
			Dorothy.realpredicate(x))

	Occupancy(x) = Dorothy.PropertySelector((model, cache) ->
			model.occupancies, Dorothy.realpredicate(x))

	SS(X...) = Dorothy.PropertySelector(Dorothy.getss,
			Dorothy.stringpredicate(X...))

	const VSite = Dorothy.CachedSelector(Name(isvsite), :system, :vsites)

	const Hydrogen = Dorothy.CachedSelector(Element("H"), :system, :hydrogens)

	const Heavy = ! (Hydrogen | VSite)

	const Water = Dorothy.CachedSelector(ResName(iswater), :system, :water)

	const Protein =
			Dorothy.CachedSelector(ResName(isprotein), :system, :protein)

	const AcidResidue = Dorothy.CachedSelector(ResName(isacidresidue),
			:system, :acidresidues)

	const BasicResidue = Dorothy.CachedSelector(ResName(isbasicresidue),
			:system, :basicresidues)

	const ChargedResidue = Dorothy.CachedSelector(ResName(ischargedresidue),
			:system, :chargedresidues)

	const PolarResidue = Dorothy.CachedSelector(ResName(ispolarresidue),
			:system, :polarresidues)

	const HydrophobicResidue =
			Dorothy.CachedSelector(ResName(ishydrophobicresidue),
			:system, :hydrophobicresidues)

	const MainChainName =
			Dorothy.CachedSelector(Name(Properties.ismainchainname),
			:system, :mainchainnames)

	const MainChain = MainChainName & Protein

	const SideChain = !MainChainName & Protein

	const BackboneName =
			Dorothy.CachedSelector(Name(Properties.isbackbonename),
			:system, :backbonenames)

	const Backbone = BackboneName & Protein

	const CAlpha = Name("CA") & Protein

	const Cα = CAlpha

	const Nter = Index(1, by=Residue, ineach=Chain) & Protein

	const Cter = Index(last, by=Residue, ineach=Chain) & Protein

	const NuclAcid = Dorothy.CachedSelector(ResName(isnuclacid),
			:system, :nuclacid)

	const Lipid = Dorothy.CachedSelector(ResName(islipid), :system, :lipid)

	const Ion = Dorothy.CachedSelector(ResName(ision), :system, :ion)

	const Helix = Dorothy.CachedSelector(SS(ishelix), :frame, :helix)

	const AlphaHelix = SS(isalphahelix)

	const Helix310 = SS(ishelix310)

	const PiHelix = SS(ispihelix)

	const Turn = SS(isturn)

	const Sheet = Dorothy.CachedSelector(SS(issheet), :frame, :sheet)

	const Strand = SS(isstrand)

	const Bridge = SS(isbridge)

	const Loop = Dorothy.CachedSelector(SS(isloop), :frame, :loop)

	const Coil = SS(iscoil)

	const Bend = SS(isbend)

end
