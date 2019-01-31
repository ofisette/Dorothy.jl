# Generic utilities

module Utils

export
		τ, emptyindices, checkindexseries, isindexseries, deleteatmap,
		splicemap, subsetsequal, substrip, wraptext, ncols, nrows, Repeated

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
	lastv = I[1]
	for i in 2:length(I)
		v = I[i]
		if ! (v > lastv)
			return false
		end
		lastv = v
	end
	true
end

function deleteatmap(n::Integer, I::AbstractVector{<:Integer})
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

function splicemap(n::Integer, range::AbstractUnitRange{<:Integer},
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

struct Repeated{T1,T2,N} <: AbstractArray{T1,N}
	val::T1
	size::T2

	Repeated{T1,T2,N}(val::T1, size::T2) where {T1,T2,N} = new(val, size)
end

Repeated(val::T1, size::T2) where {T1,T2} =
		Repeated{T1,T2,length(size)}(val, size)

Repeated(val::T1, n::Integer) where {T1} = Repeated(val, (n,))

Base.size(A::Repeated) = A.size

Base.getindex(A::Repeated, i::Int) = A.val

Base.IndexStyle(::Type{<:Repeated}) = IndexLinear()

end # module
