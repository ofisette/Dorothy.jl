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

Base.read(::DorothyIO, ::MIME"text/x-ndx", io::IO, args...; kwargs...) =
        readndx!(io, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"text/x-ndx", io::IO,
        groups::AbstractDict{<:AbstractString,<:Integer}, args...; kwargs...) =
        writendx(io, groups, args...; kwargs...)

function readndx(io::IO)
    groups = Dict{String,Vector{Int}}()
    for line in eachline(io)
        if startswith(line, "[")
            name = strip(line[2:end-1])
            indices = Int[]
            groups[title] = indices
        else
            for token in split(line)
                push!(indices, parse(Int, token))
            end
        end
    end
    groups
end

function writendx(io::IO, groups::AbstractDict{<:AbstractString,<:Integer})
    for (i, (title, indices)) in enumerate(pairs(groups))
        if i != 0
            print(io, "\n")
        end
        print(io, "[ $(title) ]\n")
        print(io, join(indices, " "))
        print(io, "\n")
    end
end

end # module
