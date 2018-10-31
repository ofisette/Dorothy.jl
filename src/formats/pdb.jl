# Protein Data Bank (PDB) molecular structure file format

module PDB

using Printf
using ..Dorothy
using ..Dorothy.Utils
using ..Dorothy.Geometry
using Formats

export register_pdb, readpdb!, writepdb

function register_pdb()
    Formats.addformat("structure/x-pdb")
    Formats.addextension("structure/x-pdb", ".pdb")
    Formats.addsignature("structure/x-pdb", "HEADER")
    Formats.addsignature("structure/x-pdb", "TITLE ")
    Formats.addsignature("structure/x-pdb", "CRYST1")
    Formats.addsignature("structure/x-pdb", "ATOM  ")
    Formats.addsignature("structure/x-pdb", "HETATM")
    Formats.addreader("structure/x-pdb", DorothyIO())
    Formats.addwriter("structure/x-pdb", DorothyIO())
end

Base.read(::DorothyIO, ::MIME"structure/x-pdb", io::IO, args...; kwargs...) =
        readpdb!(io, MolecularModel(), args...; kwargs...)

Base.read!(::DorothyIO, ::MIME"structure/x-pdb", model::MolecularModel, args...;
        kwargs...) = readpdb!(io, model, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"structure/x-pdb", io::IO,
        model::ParticleCollection, args...; kwargs...) =
        writepdb(io, model, args...; kwargs...)

pdb_particle_name(resname::AbstractString) =
        length(resname) > 3 ? resname : @sprintf(" %-3s", resname)

function readpdb!(io::IO, model::MolecularModel;
        lnoffset::AbstractArray{<:Integer,0} = zeros(Int))
    ids = get!(model, :ids, undef)
    names = get!(model, :names, undef)
    resids = get!(model, :resids, undef)
    resnames = get!(model, :resnames, undef)
    chainids = get!(model, :chainids, undef)
    R = get!(model, :R, undef)
    occupancies = get!(model, :occupancies, undef)
    bfactors = get!(model, :bfactors, undef)
    elements = get!(model, :elements, undef)
    readpdb!(io, model, ids, names, resids, resnames, chainids, R, occupancies,
            bfactors, elements, lnoffset)
end

function readpdb!(io::IO, model::MolecularModel, ids::AbstractVector{<:Integer},
        names::AbstractVector{<:AbstractString},
        resids::AbstractVector{<:Integer},
        resnames::AbstractVector{<:AbstractString},
        chainids::AbstractVector{<:AbstractString},
        R::AbstractMatrix{<:Real}, occupancies::AbstractVector{<:Real},
        bfactors::AbstractVector{<:Real},
        elements::AbstractVector{<:AbstractString},
        lnoffset::AbstractArray{<:Integer,0})

    @debug "reading PDB molecular structure"
    titlelines = String[]
    cell = nothing
    atomrecords = Tuple{Int,String}[]
    ln = lnoffset[]
    while !eof(io)
        ln += 1
        line = readline(io)
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            @debug "reading particle record" ln
            push!(atomrecords, (ln, line))
        elseif startswith(line, "TITLE")
            @debug "reading title" ln
            push!(titlelines, strip(line[11:end]))
        elseif startswith(line, "CRYST1")
            @debug "reading cell" ln
            if cell != nothing
                error("multiple cell definitions")
            else
                tokens = split(line[7:end])
                lengths = [parse(Float64, token) for token in tokens[1:3]]
                angles = [parse(Float64, token) for token in tokens[4:6]]
                cell = pbccell(lengths, deg2rad.(angles))
            end
        elseif startswith(line, "END")
            @debug "stopping" ln
            break
        else
            @debug "ignoring unknown record" ln
        end
    end
    lnoffset[] = ln

    @debug "preparing molecular model"
    nparticles = length(atomrecords)
    nparticles > 0 || error("no particle record found")
    resize!(model, nparticles)
    if length(titlelines) > 0
        model.header.title = join(titlelines, " ")
    end
    if cell != nothing
        model.header.cell = cell
    end

    @debug "parsing particle records" nparticles
    prevline = repeat("\n" , 78)
    for (i, (ln, line)) in enumerate(atomrecords)
        @debug "parsing particle record" i ln
        @debug "parsing id"
        ids[i] = parse(Int, SubString(line, 7, 11))
        @debug "parsing name"
        names[i] = substrip(line, 13, 16)
        @debug "parsing resname"
        if subsetsequal(line, prevline, 18, 21)
            resnames[i] = resnames[i-1]
        else
            resnames[i] = substrip(line, 18, 21)
        end
        @debug "parsing chainid"
        if subsetsequal(line, prevline, 22, 22)
            chainids[i] = chainids[i-1]
        else
            chainids[i] = substrip(line, 22, 22)
        end
        @debug "parsing resid"
        resids[i] = parse(Int, SubString(line, 23, 26))
        @debug "parsing Rx"
        R[1,i] = parse(Float64, SubString(line, 31, 38))
        @debug "parsing Ry"
        R[2,i] = parse(Float64, SubString(line, 39, 46))
        @debug "parsing Rz"
        R[3,i] = parse(Float64, SubString(line, 47, 54))
        @debug "parsing occupancy"
        occupancies[i] = parse(Float64, SubString(line, 55, 60))
        @debug "parsing B"
        bfactors[i] = parse(Float64, SubString(line, 61, 66))
        @debug "parsing element"
        if subsetsequal(line, prevline, 77, 78)
            elements[i] = elements[i-1]
        else
            elements[i] = substrip(line, 77, 78)
        end
        prevline = line
    end
    unwrapids!(ids, 99999)
    unwrapnames!(names)
    unwrapids!(resids, 9999)

    model
end

function writepdb(io::IO, model::ParticleCollection;
        modelid::Union{Integer,Nothing} = nothing, endrecord::Bool = true,
        lnoffset::AbstractArray{<:Integer,0} = zeros(Int))
    nparticles = length(model)
    writepdb(io, nparticles, get(model.header, :title, nothing),
            get(model.header, :cell, nothing), model.ids, model.names,
            model.resids, model.resnames,
            get(model, :chainids, ScalarArray("", nparticles)), model.R,
            get(model, :occupancies, ScalarArray(1.0, nparticles)),
            get(model, :bfactors, ScalarArray(0.0, nparticles)),
            get(model, :elements, ScalarArray("", nparticles)), modelid,
            endrecord, lnoffset)
end

function writepdb(io::IO, nparticles::Integer,
        title::Union{AbstractString,Nothing},
        cell::Union{AbstractMatrix{<:Real},Nothing},
        ids::AbstractVector{<:Integer}, names::AbstractVector{<:AbstractString},
        resids::AbstractVector{<:Integer},
        resnames::AbstractVector{<:AbstractString},
        chainids::AbstractVector{<:AbstractString}, R::AbstractMatrix{<:Real},
        occupancies::AbstractVector{<:Real}, bfactors::AbstractVector{<:Real},
        elements::AbstractVector{<:AbstractString},
        modelid::Union{Integer,Nothing}, endrecord::Bool,
        lnoffset::AbstractArray{<:Integer,0})
    @debug "writing PDB molecular structure"

    @debug "checking model properties"
    nparticles > 0 || error("number of particles must be strictly positive")
    i = findfirst(x -> length(x) > 4, names)
    i == nothing || error("particle $(i): name exceeds 4 chars")
    i = findfirst(x -> length(x) > 4, resnames)
    i == nothing || error("particle $(i): resname exceeds 4 chars")
    i = findfirst(x -> length(x) > 1, chainids)
    i == nothing || error("particle $(i): chainid exceeds 1 char")
    i = findfirst(x -> length(x) > 2, elements)
    i == nothing || error("particle $(i): element exceeds 2 chars")

    ln = lnoffset[]

    if title != nothing
        ln += 1
        @debug "writing title" ln
        for (i, line) in enumerate(wraptext(title, 70))
            if i == 1
                @printf(io, "TITLE     %s\n", line)
            else
                ln += 1
                @debug "continuing title" ln
                @printf(io, "TITLE %4d %s\n", i, line)
            end
        end
    end

    if cell != nothing
        ln += 1
        @debug "writing cell" ln
        lengths, angles = pbcbox(cell)
        @printf(io, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
                lengths..., rad2deg.(angles)...)
    end

    if modelid != nothing
        ln += 1
        @debug "start of model" ln
        @printf(io, "MODEL     %4d\n", modelid)
    end

    @debug "writing particle records" nparticles
    for i = 1:nparticles
        ln += 1
        @debug "writing particle record" i ln
        id = wrapid(ids[i], 99999)
        name = pdb_particle_name(names[i])
        resid = wrapid(resids[i], 9999)
        x, y, z = R[1,i], R[2,i], R[3,i]
        element_i = (elements[i] == "H" ? "" : elements[i])
        @printf(io, "ATOM  %5d %4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
                id, name, resnames[i], chainids[i], resid, x, y, z,
                occupancies[i], bfactors[i], elements[i])
    end

    if modelid != nothing
        ln += 1
        @debug "end of model" ln
        write(io, "ENDMDL\n")
    end

    if endrecord
        ln += 1
        @debug "end of PDB" ln
        write(io, "END\n")
    end

    lnoffset[] = ln
end

end # module
