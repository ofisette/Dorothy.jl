# Generic utilities

module Utils

using .Iterators: drop

export
		τ, emptyindices, checkindexseries, isindexseries, deleteat_map,
		splice_map, readuntil!, subsetsequal, substrip, wraptext, ncols,
		eachcol, nrows, eachrow, VectorBasedMatrix, FixedArray, ReadonlyArray,
		ScalarArray, SingleScalarVector, SelfVector, RangeVector

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

function readuntil!(src::IO, dest::AbstractArray{UInt8}, delim::UInt8)
	empty!(dest)
	while ! eof(src)
		byte = read(src, UInt8)
		if byte == delim
			break
		else
			push!(dest, byte)
		end
	end
	dest
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

struct MatrixColumns{T1<:AbstractMatrix,T2} <: AbstractVector{T2}
	M::T1

	MatrixColumns{T1,T2}(M::T1) where {T1<:AbstractMatrix,T2} = new(M)
end

MatrixColumns(M::T) where {T<:AbstractMatrix} = MatrixColumns{T,eltype(T)}(M)

Base.size(cols::MatrixColumns) = ncols(cols.M)

Base.getindex(cols::MatrixColumns, i::Int) = cols.M[:,i]

Base.setindex!(cols::MatrixColumns, v, i::Int) = cols.M[:,i] = v

Base.IndexStyle(::Type{<:MatrixColumns}) = IndexLinear()

eachcol(M::AbstractMatrix) = MatrixColumns(M)

nrows(M::AbstractMatrix) = size(M, 1)

struct MatrixRows{T1<:AbstractMatrix,T2} <: AbstractVector{T2}
	M::T1

	MatrixRows{T1,T2}(M::T1) where {T1<:AbstractMatrix,T2} = new(M)
end

MatrixRows(M::T) where {T<:AbstractMatrix} = MatrixRows{T,eltype(T)}(M)

Base.size(rows::MatrixRows) = nrows(rows.M)

Base.getindex(rows::MatrixRows, i::Int) = rows.M[i,:]

Base.setindex!(rows::MatrixRows, v, i::Int) = rows.M[i,:] = v

Base.IndexStyle(::Type{<:MatrixRows}) = IndexLinear()

eachrow(M::AbstractMatrix) = MatrixRows(M)

struct VectorBasedMatrix{T} <: AbstractMatrix{T}
	V::Vector{T}
	nrows::Int

	function VectorBasedMatrix{T}(V::Vector{T}, nrows::Integer) where {T}
		@boundscheck length(V) % nrows == 0 ||
				error("V: length must be a multiple of nrows")
		new(V, nrows)
	end
end

VectorBasedMatrix(V::Vector{T}, nrows::Integer) where {T} =
		VectorBasedMatrix{T}(V, nrows)

Base.size(M::VectorBasedMatrix) = (M.nrows, length(M.V) ÷ M.nrows)

function Base.getindex(M::VectorBasedMatrix, i::Int)
	@boundscheck checkbounds(M, i)
	@inbounds M.V[i]
end

function Base.setindex!(M::VectorBasedMatrix, v, i::Int)
	@boundscheck checkbounds(M, i)
	@inbounds M.V[i] = v
end

Base.IndexStyle(::Type{<:VectorBasedMatrix}) = IndexLinear()

function Base.empty!(M::VectorBasedMatrix)
	empty!(M.V)
	M
end

function colstoindices(cols::AbstractVector{<:Integer}, nrows::Integer)
	n = length(cols) * nrows
	I = similar(cols, length(cols) * nrows)
	i = 1
	for col in cols
		offset = (col - 1) * nrows
		for row in 1:nrows
			I[i] = offset + row
			i += 1
		end
	end
	I
end

function colstoindices(cols::UnitRange{<:Integer}, nrows::Integer)
	i = (first(cols) - 1) * nrows + 1
	j = last(cols) * nrows
	i:j
end

function Base.deleteat!(M::VectorBasedMatrix, I::AbstractVector{<:Integer})
	@boundscheck checkbounds(M, I)
	@inbounds deleteat!(M.V, colstoindices(I))
	M
end

function Base.splice!(M::VectorBasedMatrix, I::UnitRange{<:Integer},
		replacement)
	@boundscheck checkbounds(M, I)
	@inbounds splice!(M.V, colstoindices(I),
			view(replacement, length(replacement)))
end

function Base.resize!(M::VectorBasedMatrix, n::Integer)
	resize!(M.V, n * M.nrows)
	M
end

Base.sizehint!(M::VectorBasedMatrix, n::Integer) = sizehint!(M.V, n * M.nrows)

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

struct ReadonlyArray{T1,T2,N} <: AbstractArray{T2,N}
	A::T1

	ReadonlyArray{T1,T2,N}(A::T1) where {T2,N,T1<:AbstractArray{T2,N}} = new(A)
end

ReadonlyArray(A::AbstractArray{T}) where {T} =
		ReadonlyArray{typeof(A),T,ndims(A)}(A)

Base.size(A::ReadonlyArray) = size(A.A)

function Base.getindex(A::ReadonlyArray, i::Int)
	@boundscheck checkbounds(A, i)
	@inbounds getindex(A.A, i)
end

Base.IndexStyle(::Type{<:ReadonlyArray}) = IndexLinear()

struct ScalarArray{T1,T2,N} <: AbstractArray{T1,N}
	val::T1
	size::T2

	ScalarArray{T1,T2,N}(val::T1, size::T2) where {T1,T2,N} = new(val, size)
end

ScalarArray(val::T1, size::T2) where {T1,T2} =
		ScalarArray{T1,T2,length(size)}(val, size)

ScalarArray(val::T1, n::Integer) where {T1} = ScalarArray(val, (n,))

Base.size(A::ScalarArray) = A.size

Base.getindex(A::ScalarArray, i::Int) = A.val

Base.setindex!(A::ScalarArray, v, i::Int) = (A.val = v)

Base.IndexStyle(::Type{<:ScalarArray}) = IndexLinear()

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
