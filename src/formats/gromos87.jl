# Gromos87 molecular structure file format

module Gromos87

using Printf
using ..Dorothy
using ..Dorothy.Utils
using Formats

export register_gromos87, readgromos87!, writegromos87

const gromos87_cellorder =
        ((1,1), (2,2), (3,3), (2,1), (3,1), (1,2), (3,2), (1,3), (2,3))

function register_gromos87()
    Formats.addformat("structure/x-gro")
    Formats.addextension("structure/x-gro", ".gro")
    Formats.addreader("structure/x-gro", DorothyIO())
    Formats.addwriter("structure/x-gro", DorothyIO())
end

Base.read(::DorothyIO, ::MIME"structure/x-gro", io::IO, args...; kwargs...) =
        readgromos87!(io, MolecularModel(), args...; kwargs...)

Base.read!(::DorothyIO, ::MIME"structure/x-gro", io::IO, model::MolecularModel,
        args...; kwargs...) = readgromos!(io, model, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"structure/x-gro", io::IO,
        model::ParticleCollection, args...; kwargs...) =
        writegromos87(io, model, args...; kwargs...)

function readgromos87!(io::IO, model::MolecularModel;
        lnoffset::AbstractArray{<:Integer,0} = zeros(Int))
    @debug "reading Gromos87 molecular structure"

    ln = lnoffset[] + 1
    @debug "reading header" ln
    line = readline(io)
    @debug "reading title"
    tokens = rsplit(line, ", t=", limit=2)
    title = tokens[1]
    if length(tokens) == 2
        @debug "reading time"
        time = parse(Float64, tokens[2])
    else
        time = nothing
    end

    ln += 1
    @debug "reading number of particles" ln
    line = readline(io)
    nparticles = parse(Int, line)
    nparticles > 0 ||
            error("$(ln): number of particles must be strictly positive")

    @debug "reading particle records" nparticles
    records = sizehint!(Tuple{Int,String}[], nparticles)
    for i = 1:nparticles
        ln += 1
        @debug "reading particle record" i ln
        push!(records, (ln, readline(io)))
    end

    ln += 1
    @debug "reading cell" ln
    line = readline(io)
    tokens = split(line)
    3 <= length(tokens) <= 9 ||
            error("l$(ln): expected 3 to 9 cell vector components")
    cell = zeros(3,3)
    for (i, token) in enumerate(tokens)
        cell[gromos87_cellorder[i]...] = parse(Float64, token) * 10
    end
    lnoffset[] = ln

    @debug "preparing molecular model"
    resize!(model, nparticles)
    model.header.title = title
    if time != nothing
        model.header.time = time
    end
    model.header.cell = cell
    ids = get!(model, :ids, undef)
    names = get!(model, :names, undef)
    resids = get!(model, :resids, undef)
    resnames = get!(model, :resnames, undef)
    R = get!(model, :R, undef)
    ln, line = records[1]
    ncols = length(line)
    if ncols == 68
        V = get!(model, :V, undef)
    elseif ncols == 44
        V = nothing
    else
        error("$(ln): unexpected number of columns in particle record 1")
    end

    readgromos87!(ids, names, resids, resnames, R, V, records)

    model
end

function readgromos87!(ids::AbstractVector{<:Integer},
        names::AbstractVector{<:AbstractString},
        resids::AbstractVector{<:Integer},
        resnames::AbstractVector{<:AbstractString}, R::AbstractMatrix{<:Real},
        V::Union{AbstractMatrix{<:Real},Nothing},
        records::AbstractVector{<:Tuple{Integer,AbstractString}})
    prevline = repeat(" " , 68)
    for (i, (ln, line)) in enumerate(records)
        @debug "parsing particle record" i ln
        @debug "parsing resid"
        resids[i] = parse(Int, SubString(line, 1, 5))
        @debug "parsing resname"
        if subsetsequal(line, prevline, 6, 10)
            resnames[i] = resnames[i-1]
        else
            resnames[i] = strip(line, 6, 10)
        end
        @debug "parsing name"
        names[i] = strip(line, 11, 15)
        @debug "parsing id"
        ids[i] = parse(Int, SubString(line, 16, 20))
        @debug "parsing Rx"
        R[1,i] = parse(Float64, SubString(line, 21, 28)) * 10.0
        @debug "parsing Ry"
        R[2,i] = parse(Float64, SubString(line, 29, 36)) * 10.0
        @debug "parsing Rz"
        R[3,i] = parse(Float64, SubString(line, 37, 44)) * 10.0
        if V != nothing
            @debug "parsing Vx"
            V[1,i] = parse(Float64, SubString(line, 45, 52)) * 10.0
            @debug "parsing Vy"
            V[2,i] = parse(Float64, SubString(line, 53, 60)) * 10.0
            @debug "parsing Vz"
            V[3,i] = parse(Float64, SubString(line, 61, 68)) * 10.0
        end
        prevline = line
    end
    unwrapids!(ids, 99999)
    unwrapids!(resids, 99999)
end

writegromos87(io::IO, model::ParticleCollection;
        lnoffset::AbstractArray{<:Integer,0} = zeros(Int)) =
        writegromos87(io, length(model), get(model.header, :title, ""),
        get(model.header, :time, nothing), get(model.header, :cell, zeros(3,3)),
        model.ids, model.names, model.resids, model.resnames, model.R,
        get(model, :V, nothing), lnoffset)

function writegromos87(io::IO, nparticles::Integer, title::AbstractString,
        time::Union{Real,Nothing}, cell::AbstractMatrix{<:Real},
        ids::AbstractVector{<:Integer}, names::AbstractVector{<:AbstractString},
        resids::AbstractVector{<:Integer},
        resnames::AbstractVector{<:AbstractString}, R::AbstractMatrix{<:Real},
        V::Union{AbstractMatrix{<:Real},Nothing},
        lnoffset::AbstractArray{<:Integer,0})
    @debug "writing Gromos87 molecular structure"

    @debug "checking model properties"
    nparticles > 0 || error("number of particles must be strictly positive")
    i = findfirst(x -> length(x) > 5, names)
    i == nothing || error("particle $(i): name exceeds 5 chars")
    i = findfirst(x -> length(x) > 5, resnames)
    i == nothing || error("particle $(i): resname exceeds 5 chars")

    ln = lnoffset[] + 1
    @debug "writing header" ln
    header = title
    if time != nothing
        header = header * ", t= $(time)"
    end
    write(io, header * "\n")

    ln += 1
    @debug "writing number of particles" ln
    write(io, "$(nparticles)\n")

    @debug "writing particle records" nparticles
    for i = 1:nparticles
        ln += 1
        @debug "writing particle record" i ln
        if V == nothing
            @printf(io, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                    wrapid(resids[i], 99999), resnames[i],
                    names[i], wrapid(ids[i], 99999),
                    R[1,i] / 10.0, R[2,i] / 10.0, R[3,i] / 10.0)
        else
            @printf(io, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                    wrapid(resids[i], 99999), resnames[i],
                    names[i], wrapid(ids[i], 99999),
                    R[1,i] / 10.0, R[2,i] / 10.0, R[3,i] / 10.0,
                    V[1,i] / 10.0, V[2,i] / 10.0, V[3,i] / 10.0)
        end
    end

    ln += 1
    @debug "writing cell" ln
    gromos87_cell = [cell[i,j] / 10.0 for (i,j) in gromos87_cellorder]
    if all(x -> isapprox(x, 0.0, atol=0.00001), gromos87_cell[4:9])
        @printf(io, "   %.5f   %.5f   %.5f\n", gromos87_cell[1:3]...)
    else
        @printf(io,
            "   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f   %.5f\n",
            gromos87_cell...)
    end
    lnoffset[] = ln
end

end # module
