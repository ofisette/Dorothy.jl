# Gromacs index files

module NDX

using ..Dorothy
using Formats

export register_ndx, readndx, writendx

function register_ndx()
	Formats.addformat("text/x-ndx")
	Formats.addextension("text/x-ndx", ".ndx")
	Formats.addreader("text/x-ndx", DorothyIO())
	Formats.addwriter("text/x-ndx", DorothyIO())
end

const AbstractIndex = AbstractDict{<:AbstractString,<:AbstractVector{<:Integer}}

Base.read(::DorothyIO, ::MIME"text/x-ndx", io::IO, args...; kwargs...) =
		readndx!(io, Dict{String,Vector{Int}}(), args...; kwargs...)

Base.read!(::DorothyIO, ::MIME"text/x-ndx", io::IO, index::AbstractIndex,
		args...; kwargs...) = readndx!(io, index, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"text/x-ndx", io::IO, index::AbstractIndex,
		args...; kwargs...) = writendx(io, index, args...; kwargs...)

function readndx!(io::IO, index::AbstractIndex)
	empty!(index)
	local I::Vector{Int}
	for line in eachline(io)
		if startswith(line, "[")
			group = strip(line[2:end-1])
			I = Int[]
			index[group] = I
		else
			for token in split(line)
				push!(I, parse(Int, token))
			end
		end
	end
	index
end

function writendx(io::IO, index::AbstractIndex)
	for (n, (group, I)) in enumerate(pairs(index))
		if n != 0
			print(io, "\n")
		end
		print(io, "[ $(group) ]\n")
		print(io, join(I, " "))
		print(io, "\n")
	end
end

end # module
