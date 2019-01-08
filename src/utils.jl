# Generic utilities

module Utils

using .Iterators: drop

export
		τ, emptyindices, checkindexseries, isindexseries, deleteat_map,
		splice_map, subsetsequal, substrip, wraptext, ncols, nrows, FixedArray,
		Repeat, SingleScalarVector, SelfVector, RangeVector

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

function splice_map(n::Integer, range::UnitRange{<:Integer},
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

Base.setindex!(A::Repeat, v, i::Int) = (A.val = v)

Base.IndexStyle(::Type{<:Repeat}) = IndexLinear()

struct SingleScalarVector{T} <: AbstractVector{T}
	val::T

	SingleScalarVector{T}(val::T) where {T} = new(val)
end

SingleScalarVector(val::T) where {T} = SingleScalarVector{T}(val)

Base.size(V::SingleScalarVector) = (1,)

Base.getindex(V::SingleScalarVector, i::Int) = V.val

Base.setindex!(V::SingleScalarVector, v, i::Int) = (V.val = v)

Base.IndexStyle(::Type{<:SingleScalarVector}) = IndexLinear()

struct SelfVector <: AbstractVector{SingleScalarVector{Int}}
	n::Int

	SelfVector(n::Integer) = new(n)
end

Base.size(V::SelfVector) = (V.n,)

Base.getindex(V::SelfVector, i::Int) = SingleScalarVector(i)

Base.IndexStyle(::Type{<:SelfVector}) = IndexLinear()

struct RangeVector{T<:Integer} <: AbstractVector{T}
	ranges::Vector{UnitRange{T}}
	n::Int

	function RangeVector{T}(ranges::Vector{UnitRange{T}}) where {T<:Integer}
		@boundscheck !isempty(ranges) || error("expected at least one range")
		new(ranges, sum(length.(ranges)))
	end
end

RangeVector(ranges::Vector{UnitRange{T}}) where {T<:Integer} =
		RangeVector{T}(ranges)

Base.size(V::RangeVector) = (V.n,)

function Base.getindex(V::RangeVector, i::Integer)
	@boundscheck checkbounds(V, i)
	offset = 0
	for r in V.ranges
		n = length(r)
		if i > n + offset
			offset += n
		else
			return r[i-offset]
		end
	end
end

Base.IndexStyle(::Type{<:RangeVector}) = IndexLinear()

end # module
