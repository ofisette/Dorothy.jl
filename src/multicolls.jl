# Type for collections of arrays with restricted choice of keys and value types

module Multicollections

using ..Dorothy.Utils
using ..Dorothy.Graphs

export
        AbstractMulticollection, Multicollection, MulticollectionView,
        MulticollectionItem, collval, itemtocollprop, colltoitemprop,
        collitemname, MulticollectionSplit, eachitem, asawhole

abstract type AbstractMulticollection end

abstract type Multicollection <: AbstractMulticollection
    #=
    mutable struct
        header<:Header
        n::Int
        D<:Dict{Symbol,Any}
    end
    =#
end

struct MulticollectionView{T1<:Multicollection,
        T2<:AbstractVector{<:Integer}}  <: AbstractMulticollection
    C::T1
    I::T2

    @inline function MulticollectionView{T1,T2}(C::T1, I::T2) where
            {T1<:Multicollection,T2<:AbstractVector{<:Integer}}
        @boundscheck checkindexseries(C, I)
        new(C, I)
    end
end

struct MulticollectionItem{T<:Multicollection}
    C::T
    i::Int

    @inline function MulticollectionItem{T}(C::T, i::Integer) where
            {T<:Multicollection}
        @boundscheck checkbounds(C, i)
        new(C, i)
    end
end

@inline MulticollectionView(C::T1, I::T2) where
        {T1<:Multicollection,T2<:AbstractVector{<:Integer}} =
        MulticollectionView{T1,T2}(C, I)

@inline function MulticollectionView(C::MulticollectionView,
        I::AbstractVector{<:Integer})
    @boundscheck checkindexseries(C, I)
    @inbounds MulticollectionView(C.C, C.I[I])
end

MulticollectionItem(C::T, i::Integer) where {T<:Multicollection} =
        MulticollectionItem{T}(C, i)

MulticollectionItem(C::MulticollectionView{T}, i::Integer) where
        {T<:Multicollection} = MulticollectionItem{T}(C.C, C.I[i])

@inline Base.getproperty(C::AbstractMulticollection, name::Symbol) =
		getproperty(C, Val(name))

@inline Base.getproperty(C::AbstractMulticollection, ::Val{name}) where {name} =
		getindex(C, name)

@inline Base.getproperty(C::Multicollection, ::Val{:header}) =
		getfield(C, :header)

@inline Base.getproperty(C::Multicollection, ::Val{:n}) = getfield(C, :n)

@inline Base.getproperty(C::Multicollection, ::Val{:D}) = getfield(C, :D)

@inline Base.getproperty(C::MulticollectionView, ::Val{:C}) = getfield(C, :C)

@inline Base.getproperty(C::MulticollectionView, ::Val{:I}) = getfield(C, :I)

@inline Base.getproperty(C::MulticollectionView, ::Val{:header}) = C.C.header

@inline Base.getproperty(item::MulticollectionItem, name::Symbol) =
		Base.getproperty(item, Val(name))

@inline Base.getproperty(item::MulticollectionItem, ::Val{name}) where {name} =
		getindex(item, name)

@inline Base.getproperty(item::MulticollectionItem, ::Val{:C}) =
        getfield(item, :C)

@inline Base.getproperty(item::MulticollectionItem, ::Val{:i}) =
        getfield(item, :i)

@inline Base.getproperty(item::MulticollectionItem, ::Val{:header}) =
        item.C.header

@inline Base.setproperty!(C::AbstractMulticollection, name::Symbol, v) =
        setproperty!(C, Val(name), v)

@inline Base.setproperty!(C::AbstractMulticollection, ::Val{name}, v) where
        {name} = setindex!(C, v, name)

@inline Base.setproperty!(C::Multicollection, ::Val{:n}, v) =
        setfield!(C, :n, v)

@inline Base.setproperty!(item::MulticollectionItem, name::Symbol, v) =
        setindex!(item, v, name)

Base.propertynames(C::Multicollection, private::Bool = false) =
        private ? (:header, :n, :D, keys(C)...) : (:header, keys(C)...)

Base.propertynames(C::MulticollectionView, private::Bool = false) =
		private ? (:header, :C, :I, keys(C.C)...) : (:header, keys(C.C)...)

function Base.propertynames(item::MulticollectionItem, private::Bool = false)
    dynprops = [colltoitemprop(item, key) for key in keys(item.C)]
	private ? (:header, :C, :i, dynprops...) : (:header, dynprops...)
end

Base.show(io::IO, C::AbstractMulticollection) =
        print(io, "$(typeof(C))($(length(C)))")

function Base.show(io::IO, ::MIME"text/plain", C::Multicollection)
	n = length(C)
	nkeys = length(keys(C))
    itemname = collitemname(C)
	print(io, "$(n)-$(itemname) $(nkeys)-key $(typeof(C))")
	if nkeys > 0
		print(io, ":")
		for key in sort!(collect(keys(C)), by = lowercase∘string)
			print(io, "\n $(key)")
		end
	end
end

function Base.show(io::IO, ::MIME"text/plain", C::MulticollectionView)
    n = length(C)
    nkeys = length(keys(C))
    itemname = collitemname(C.C)
	print(io, "$(n)-$(itemname) $(nkeys)-key $(typeof(C.C)) view")
    if nkeys > 0
        print(io, ":")
        for key in sort!(collect(keys(C)), by = lowercase∘string)
            print(io, "\n $(key)")
        end
    end
end

Base.show(io::IO, item::MulticollectionItem) =
        print(io, "$(typeof(item.C))[$(item.i)]")

function Base.show(io::IO, ::MIME"text/plain", item::MulticollectionItem)
    n = length(item.C)
    itemkeys = keys(item)
    nkeys = length(itemkeys)
    itemname = collitemname(item.C)
	print(io, "$(nkeys)-key $(typeof(item.C)) $(itemname)")
    if nkeys > 0
		print(io, ":")
		for key in sort!(itemkeys, by = lowercase∘string)
			print(io, "\n $(key)")
		end
	end
end

Base.length(C::Multicollection) = C.n

Base.length(C::MulticollectionView) = length(C.I)

Base.firstindex(C::AbstractMulticollection) = 1

Base.lastindex(C::AbstractMulticollection) = length(C)

Base.eachindex(C::AbstractMulticollection) = Base.OneTo(length(C))

Base.checkbounds(::Type{Bool}, C::AbstractMulticollection, I) =
		checkbounds(Bool, 1:length(C), I)

Base.checkbounds(C::AbstractMulticollection, I) = checkbounds(1:length(C), I)

function Base.iterate(C::AbstractMulticollection, i::Int = 1)
	if i > length(C)
		nothing
	else
		MulticollectionItem(C, i), i + 1
	end
end

Base.eltype(::Type{T}) where {T<:AbstractMulticollection} =
        MulticollectionItem{T}

Base.similar(C::AbstractMulticollection, t::Type, n::Integer = length(C)) =
        similar(eachindex(C), t, n)

Base.parent(C::Multicollection) = C

Base.parent(C::MulticollectionView) = C.C

Base.parent(item::MulticollectionItem) = item.C

Base.parentindices(C::Multicollection) = eachindex(C)

Base.parentindices(C::MulticollectionView) = C.I

Base.parentindices(item::MulticollectionItem) = [item.i]

Base.view(C::Multicollection, I::AbstractVector{<:Integer}) =
		MulticollectionView(C, copy(I))

Base.view(C::MulticollectionView, I::AbstractVector{<:Integer}) =
        MulticollectionView(C, I)

Base.view(C::AbstractMulticollection, i::Integer) = MulticollectionView(C, [i])

@inline Base.getindex(C::AbstractMulticollection, i::Integer) =
        MulticollectionItem(C, i)

@inline Base.getindex(C::Multicollection, key::Symbol) = C.D[key]

@inline Base.getindex(C::MulticollectionView, key::Symbol) =
        collview(C.C[key], C.I)

@inline Base.getindex(item::MulticollectionItem, key::Symbol) =
        collscalar(item.C[itemtocollprop(item, key)], item.i)

@inline function Base.setindex!(C::Multicollection, v, key::Symbol)
    if haskey(C, key)
        C[key] .= v
    else
        get!(C, key, v)
    end
end

@inline Base.setindex!(C::MulticollectionView, v, key::Symbol) = (C[key] .= v)

@inline Base.setindex!(item::MulticollectionItem, v, key::Symbol) =
        (collview(item.C[itemtocollprop(item, key)], [item.i]) .= v)

Base.haskey(C::Multicollection, key::Symbol) = haskey(C.D, key)

Base.haskey(C::MulticollectionView, key::Symbol) = haskey(C.C, key)

Base.haskey(item::MulticollectionItem, key::Symbol) =
        haskey(item.C, itemtocollprop(item, key))

Base.getkey(C::Multicollection, key::Symbol, default) =
		getkey(C.D, key, default)

Base.getkey(C::MulticollectionView, key::Symbol, default) =
		getkey(C.C, key, default)

Base.getkey(item::MulticollectionItem, key::Symbol, default) =
        haskey(item, key) ? key : default

Base.get(f, C::AbstractMulticollection, key::Symbol) =
        haskey(C, key) ? getindex(C, key) : f()

Base.get(C::AbstractMulticollection, key::Symbol, default) =
        (get(C, key) do; default; end)

Base.get(f, item::MulticollectionItem, key::Symbol) =
        haskey(item, key) ? getindex(item, key) : f()

Base.get(item::MulticollectionItem, key::Symbol, default) =
        (get(item, key) do; default; end)

function Base.get!(f, C::Multicollection, key::Symbol)
	get!(C.D, key) do
        collval(C, Val(key), length(C), f())
	end
end

Base.get!(C::Multicollection, key::Symbol, default) =
		(get!(C, key) do; default; end)

function Base.delete!(C::Multicollection, key::Symbol)
    delete!(C.D, key)
    C
end

Base.keys(C::Multicollection) = keys(C.D)

Base.keys(C::MulticollectionView) = keys(C.C)

Base.keys(item::MulticollectionItem) =
        [colltoitemprop(item, key) for key in keys(item.C)]

Base.values(C::AbstractMulticollection) = [C[key] for key in keys(C)]

Base.pairs(C::AbstractMulticollection) = [(key, C[key]) for key in keys(C)]

Base.:(==)(C1::Multicollection, C2::Multicollection) = (C1 === C2)

Base.:(==)(C1::MulticollectionView, C2::MulticollectionView) =
		(C1.C == C2.C && C1.I == C2.I)

Base.:(==)(item1::MulticollectionItem, item2::MulticollectionItem) =
        (item1.C == item2.C && item1.i == item2.i)

function Base.resize!(C::Multicollection, n::Integer)
    @boundscheck n >= 0 || error("expected positive length")
    for val in values(C.D)
        resize!(collunwrap(val), n)
    end
    C.n = n
    C
end

function Base.empty!(C::Multicollection)
    for val in values(C.D)
        empty!(collunwrap(val))
    end
    C.n = 0
    C
end

Base.deleteat!(C::Multicollection, i::Integer) = deleteat!(C, [i])

@inline function Base.deleteat!(C::Multicollection, I::AbstractArray{<:Integer})
    @boundscheck checkindexseries(C, I)
    for val in values(C.D)
        deleteat!(collunwrap(val), I)
    end
    C.n -= length(I)
    C
end

Base.splice!(C::Multicollection, i::Integer) = deleteat!(C, i)

Base.splice!(C::Multicollection, i::Integer, replacement::MulticollectionItem) =
        splice!(C, i:i, view(parent(replacement), parentindices(replacement)))

Base.splice!(C::Multicollection, i::Integer,
        replacement::AbstractMulticollection) = splice!(C, i:i, replacement)

Base.splice!(C::Multicollection, range::UnitRange{<:Integer}) =
        deleteat!(C, range)

Base.splice!(C::Multicollection, range::UnitRange{<:Integer},
        replacement::MulticollectionItem) =
        splice!(C, range, view(parent(replacement), parentindices(replacement)))

function Base.splice!(C::Multicollection, range::UnitRange{<:Integer},
        replacement::AbstractMulticollection)
    @boundscheck begin
        checkbounds(C, range)
        parent(replacement) !== C ||
                error("replacement must belong to a different $(typeof(C))")
        for key in keys(C)
            if ! haskey(replacement, key)
                error("missing $(key) key in replacement")
            end
        end
    end
    nreplaced = length(range)
    nreplacement = length(replacement)
    spliced = similar(C, nreplaced)
    if nreplacement == 0
		@inbounds deleteat!(C, range)
	else
        replaced = view(C, range)
        for (key, val) in pairs(C.D)
            get!(spliced, key, replaced[key])
            splice!(collunwrap(val), range, replacement[key])
        end
    end
    C.n += nreplacement - nreplaced
    spliced
end

Base.append!(C::Multicollection, src::MulticollectionItem) =
        append!(C, view(parent(src), parentindices(src)))

function Base.append!(C::Multicollection, src::AbstractMulticollection)
    n = length(C)
    splice!(C, n+1:n, src)
    C
end

Base.prepend!(C::Multicollection, src::MulticollectionItem) =
        prepend!(C, view(parent(src), parentindices(src)))

function Base.prepend!(C::Multicollection, src::AbstractMulticollection)
    splice!(C, 1:0, src)
    C
end

Base.insert!(C::Multicollection, i::Integer, src::MulticollectionItem) =
        insert!(C, i, view(parent(src), parentindices(src)))

function Base.insert!(C::Multicollection, i::Integer,
        src::AbstractMulticollection)
    splice!(C, i:i-1, src)
    C
end

collview(X, I::AbstractArray{<:Integer}) = view(X, I)

collview(M::AbstractMatrix, I::AbstractArray{<:Integer}) = view(M, :, I)

collscalar(X, i::Integer) = X[i]

collscalar(M::AbstractMatrix, i::Integer) = view(M, :, i)

collunwrap(x) = x

collunwrap(A::FixedArray) = A.A

collunwrap(G::FixedGraph) = G.G

function collval end

@inline itemtocollprop(item::MulticollectionItem, name::Symbol) =
        itemtocollprop(item, Val(name))

@inline itemtocollprop(::MulticollectionItem, ::Val{name}) where name = name

@inline colltoitemprop(item::MulticollectionItem, name::Symbol) =
        colltoitemprop(item, Val(name))

collitemname(::Multicollection) = "item"

@inline colltoitemprop(::MulticollectionItem, ::Val{name}) where name = name

struct MulticollectionSplit{T1<:AbstractMulticollection,
        T2<:AbstractVector{<:AbstractVector{<:Integer}},
        T3<:AbstractVector{<:Integer}}
	C::T1
	itemindices::T2
	groupindices::T3

	MulticollectionSplit{T1,T2,T3}(C::T1, itemindices::T2,
            groupindices::T3) where {T1<:AbstractMulticollection,
            T2<:AbstractVector{<:AbstractVector{<:Integer}},
            T3<:AbstractVector{<:Integer}} = new(C, itemindices, groupindices)
end

MulticollectionSplit(C::T1, itemindices::T2, groupindices::T3) where
        {T1<:AbstractMulticollection,
        T2<:AbstractVector{<:AbstractVector{<:Integer}},
        T3<:AbstractVector{<:Integer}} =
		MulticollectionSplit{T1,T2,T3}(C, itemindices, groupindices)

Base.show(io::IO, C::MulticollectionSplit) =
        print(io, "$(typeof(C))($(length(C)))")

function Base.show(io::IO, ::MIME"text/plain", C::MulticollectionSplit)
	print(io, "$(length(C))-group split of ")
    show(io, MIME"text/plain"(), C.C)
end

Base.length(split::MulticollectionSplit) = length(split.itemindices)

Base.firstindex(split::MulticollectionSplit) = 1

Base.lastindex(split::MulticollectionSplit) = length(split)

Base.eachindex(split::MulticollectionSplit) = Base.OneTo(length(split))

Base.checkbounds(::Type{Bool}, split::MulticollectionSplit, I) =
		checkbounds(Bool, 1:length(split), I)

Base.checkbounds(split::MulticollectionSplit, I) =
		checkbounds(1:length(split), I)

function Base.iterate(split::MulticollectionSplit, i::Integer = 1)
	if i > length(split.itemindices)
		nothing
	else
		MulticollectionView(split.C, split.itemindices[i]), i+1
	end
end

Base.eltype(::Type{<:MulticollectionSplit{T1,T2,T3}}) where {T1,T2,T3} =
        eltype(T2)

function Base.getindex(split::MulticollectionSplit, i::Integer)
	@boundscheck firstindex(split) <= i <= lastindex(split) ||
			throw(BoundsError(split, i))
	MulticollectionView(split.C, split.itemindices[i])
end

eachitem(C::AbstractMulticollection) =
        MulticollectionSplit(C, SelfVector(length(C)), eachindex(C))

asawhole(C::AbstractMulticollection) =
        MulticollectionSplit(C, [eachindex(C)], ScalarVector(1, length(C)))

end # module
