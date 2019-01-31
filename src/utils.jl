# Generic utilities

module Utils

using .Iterators: drop

export
		τ, emptyindices, checkindexseries, isindexseries, deleteat_map,
		splice_map, subsetsequal, substrip, wraptext, ncols, nrows, FixedArray,
		Repeat, OffsetVector, RangeIdentity

const τ = 2π

const emptyindices = Int[]

function checkindexseries(collection, I::AbstractVector{<:Integer})
	if ! isindexseries(I)
		error("expected sorted unique indices")
	end
	checkbounds(collection, I)
end

function isindexseries(I::AbstractVector{<:Integer})
	if isempty(I)
		return true
	end
	vprev = first(I)
	for v in drop(I, 1)
		if ! (v > vprev)
			return false
		end
		vprev = v
	end
	true
end

function deleteat_map(n::Integer, I::AbstractVector{<:Integer})
	map = collect(1:n)
	for i in I
		map[i] = 0
	end
	offset = 0
	for i in 1:n
		if map[i] == 0
			offset -= 1
		else
			map[i] += offset
		end
	end
	map
end

function splice_map(n::Integer, range::AbstractUnitRange{<:Integer},
		nreplacement::Integer)
	map = collect(1:n)
	for i in range
		map[i] = 0
	end
	offset = nreplacement - length(range)
	for i in last(range)+1:n
		map[i] += offset
	end
	map
end

function substrip(str::AbstractString, i0::Integer, i1::Integer)
	while i0 <= i1 && str[i0] == ' '
		i0 += 1
	end
	while i1 > i0 && str[i1] == ' '
		i1 -= 1
	end
	str[i0:i1]
end

function subsetsequal(A, B, i0, i1)
	for i = i0:i1
		if A[i] != B[i]
			return false
		end
	end
	true
end

function wraptext(text::AbstractString, maxcols::Integer = 80)
	maxcols > 0 || error("maxcols: expected strictly positive value")
	text = strip(text)
	if length(text) <= maxcols
		return [text]
	end
	tokens = split(text)
	lines = Vector{String}[]
	linecols = maxcols
	for token in tokens
		tokencols = length(token)
		if linecols + tokencols + 1 > maxcols
			push!(lines, [token])
			linecols = tokencols
		else
			push!(lines[end], token)
			linecols += tokencols + 1
		end
	end
	[join(line, " ") for line in lines]
end

ncols(M::AbstractMatrix) = size(M, 2)

nrows(M::AbstractMatrix) = size(M, 1)

struct FixedArray{T1,T2,N} <: AbstractArray{T2,N}
	A::T1

	FixedArray{T1,T2,N}(A::T1) where {T2,N,T1<:AbstractArray{T2,N}} = new(A)
end

FixedArray(A::AbstractArray{T}) where {T} = FixedArray{typeof(A),T,ndims(A)}(A)

Base.size(A::FixedArray) = size(A.A)

function Base.getindex(A::FixedArray, i::Int)
	@boundscheck checkbounds(A, i)
	@inbounds A.A[i]
end

function Base.setindex!(A::FixedArray, v, i::Int)
	@boundscheck checkbounds(A, i)
	@inbounds A.A[i] = v
end

Base.IndexStyle(::Type{<:FixedArray}) = IndexLinear()

struct Repeat{T1,T2,N} <: AbstractArray{T1,N}
	val::T1
	size::T2

	Repeat{T1,T2,N}(val::T1, size::T2) where {T1,T2,N} = new(val, size)
end

Repeat(val::T1, size::T2) where {T1,T2} =
		Repeat{T1,T2,length(size)}(val, size)

Repeat(val::T1, n::Integer) where {T1} = Repeat(val, (n,))

Base.size(A::Repeat) = A.size

Base.getindex(A::Repeat, i::Int) = A.val

Base.IndexStyle(::Type{<:Repeat}) = IndexLinear()

struct OffsetVector{T1,T2<:AbstractVector{<:T1}} <: AbstractVector{T1}
	V::T2
	d::Int

	OffsetVector{T1,T2}(V::T2, d::Integer) where {T1,T2<:AbstractVector{<:T1}} =
			new(V, d)
end

OffsetVector(V::T, d::Integer) where {T<:AbstractVector} =
		OffsetVector{eltype(T),T}(V,d)

Base.size(V::OffsetVector) = size(V.V)

Base.getindex(V::OffsetVector, i::Int) = getindex(V.V, i)

Base.setindex!(V::OffsetVector{T1,T2}, v::T1, i::Int) where
		{T1,T2<:AbstractVector{<:T1}} = setindex!(V, v, i)

Base.IndexStyle(::Type{<:OffsetVector}) = IndexLinear()

struct RangeIdentity <: AbstractVector{UnitRange{Int}}
	n::Int

	RangeIdentity(n::Integer) = new(n)
end

Base.size(V::RangeIdentity) = (V.n,)

Base.getindex(V::RangeIdentity, i::Integer) = i:i

Base.IndexStyle(::Type{<:RangeIdentity}) = IndexLinear()

end # module
