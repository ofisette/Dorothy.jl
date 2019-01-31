# Header type with restricted choice of keys and value types

module Headers

export Header, headerval

abstract type Header
	#=
	struct
		D<:Dict{Symbol,Any}
	end
	=#
end

@inline Base.getproperty(header::Header, name::Symbol) =
		getproperty(header, Val(name))

@inline Base.getproperty(header::Header, ::Val{:D}) = getfield(header, :D)

@inline Base.getproperty(header::Header, ::Val{name}) where {name} =
		getindex(header, name)

@inline Base.setproperty!(header::Header, name::Symbol, v) =
		setindex!(header, v, name)

Base.propertynames(header::Header, private::Bool = false) =
		private ? (:D, keys(header)...) : (keys(header)...,)

Base.show(io::IO, header::Header) = print(io, "$(typeof(header))()")

function Base.show(io::IO, ::MIME"text/plain", header::Header)
	print(io, "$(length(header))-entry $(typeof(header))")
	if length(header) > 0
		print(io, ":")
		for key in keys(header)
			print(io, "\n $(key)")
		end
	end
end

Base.getindex(header::Header, key::Symbol) = header.D[key]

Base.setindex!(header::Header, v, key::Symbol) =
		(header.D[key] = headerval(header, Val(key), v))

Base.haskey(header::Header, key::Symbol) = haskey(header.D, key)

Base.getkey(header::Header, key::Symbol, default) =
		getkey(header.D, key, default)

Base.get(header::Header, key::Symbol, default) = get(header.D, key, default)

Base.get(f, header::Header, key::Symbol) = get(f, header.D, key)

function Base.get!(f, header::Header, key::Symbol)
	get!(header.D, key) do
		headerval(header, Val(key), f())
	end
end

Base.get!(header::Header, key::Symbol, default) =
		(get!(header, key) do; default; end)

function Base.delete!(header::Header, key::Symbol)
	delete!(header.D, key)
	header
end

function Base.empty!(header::Header)
	empty!(header.D)
	header
end

Base.iterate(header::Header) = iterate(header.D)

Base.iterate(header::Header, state) = iterate(header.D, state)

Base.IteratorSize(header::Header) = Base.IteratorSize(header.D)

Base.IteratorEltype(header::Header) = Base.IteratorEltype(header.D)

Base.eltype(header::Header) = Pair{Symbol,Any}

Base.length(header::Header) = length(header.D)

Base.keys(header::Header) = keys(header.D)

Base.values(header::Header) = values(header.D)

Base.pairs(header::Header) = pairs(header.D)

Base.:(==)(header1::Header, header2::Header) = (header1 === header2)

function Base.merge!(header::Header, src)
	for (key, value) in pairs(src)
		header[key] = value
	end
	header
end

function headerval end

end # module
