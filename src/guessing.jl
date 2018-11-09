# Algorithms to assign particle properties

function guesselement(name::AbstractString, resname::AbstractString)
    if ismonatomicion(resname)
        titlecase(resname)
    elseif isempty(name)
        error("could not guess element from name or resname")
    else
        name[1:1]
    end
end

guesselements!(model::ParticleCollection) =
        guesselements!(get!(model, :elements, undef), model)

guesselements!(elements::AbstractVector{<:AbstractString},
        model::ParticleCollection) = guesselements!(elements, model.names,
        model.resnames)

function guesselements!(elements::AbstractVector{<:AbstractString},
        names::AbstractVector{<:AbstractString},
        resnames::AbstractVector{<:AbstractString})
    for i in eachindex(elements)
        elements[i] = guesselement(names[i], resnames[i])
    end
    elements
end

function guessmass(name::AbstractString, element::AbstractString)
    get(standard_atomic_weight, element) do
        if isvsite(name)
            0.0
        else
            error("could not guess mass from name or element")
        end
    end
end

guessmasses!(model::ParticleCollection) =
        guessmasses!(get!(model, :masses, undef), model)

guessmasses!(masses::AbstractVector{<:Real}, model::ParticleCollection) =
        guessmasses!(masses, model.names, model.elements)

function guessmasses!(masses::AbstractVector{<:Real},
        names::AbstractVector{<:AbstractString},
		elements::AbstractVector{<:AbstractString})
	for i in eachindex(masses)
        masses[i] = guessmass(names[i], elements[i])
    end
    masses
end

@inline guessss!(model::ParticleCollection; kwargs...) =
        guessss!(get!(model, :SS, ""), model; kwargs...)

@inline guessss!(ss::AbstractVector{<:AbstractString},
        model::ParticleCollection; alg::Symbol = :stride, kwargs...) =
        guessss!(Val(alg), ss, model; kwargs...)

@inline guessss!(::Val{alg}, ss::AbstractVector{<:AbstractString},
        model::ParticleCollection; kwargs...) where {alg} =
        error("unknown secondary structure algorithm: $(alg)")

@inline guessss!(::Val{:stride}, ss::AbstractVector{<:AbstractString},
        model::ParticleCollection; kwargs...) = stride!(ss, model; kwargs...)

function stride!(ss::AbstractVector{<:AbstractString},
        model::ParticleCollection; kwargs...)
    I = findall(map(Protein, model))
    protein = view(model, I)
    stride!(view(ss, I), protein.resids, protein.resnames, protein; kwargs...)
    ss
end

function stride!(ss::AbstractVector{<:AbstractString},
        resids::AbstractVector{<:Integer},
        resnames::AbstractVector{<:AbstractString}, protein::ParticleCollection;
        stridepath::AbstractString = "stride")
    itemindices = eachresidue(protein).itemindices
    striderecords = Tuple{Int,String,String}[]
    mktemp() do pdbpath, io
        write(specify(io, "structure/x-pdb"), protein)
        close(io)
        prevss = ""
        for line in readlines(open(`$(stridepath) $(pdbpath)`))
            if startswith(line, "ASG")
                tokens = split(line)
                resname = tokens[2]
                resid = parse(Int, tokens[4])
                if tokens[6] == prevss
                    ssi = prevss
                else
                    ssi = tokens[6]
                end
                push!(striderecords, (resid, resname, ssi))
                prevss = ssi
            end
        end
    end
    nstride = length(striderecords)
    nres = length(itemindices)
    nstride == nres ||
            error("cannot assign $(nstride) stride records to $(nres) residues")
    for (I, (resid, resname, ssi)) in zip(itemindices, striderecords)
        resids[I[1]] == resid ||
                error("stride resid mismatch at $(resid[I[1]])")
        resnames[I[1]] == resname ||
                error("stride resname mismatch at resid $(resname[I[1]])")
        for i in I
            ss[i] = ssi
        end
    end
    ss
end

guesstopology!(model::ParticleCollection; kwargs...) =
        guesstopology!(get!(model, :topology, []), model; kwargs...)

function guesstopology!(topology::AbstractGraph, model::ParticleCollection;
        tol::Real = 0.1)
    @boundscheck tol > 0 || error("expected strictly positive tolerance")
    cell = get(model.header, :cell) do
        lengths = dim(extent(model)) .+
                (1.0+tol)^2 * 2.0 * maximum(values(covalent_radius))
        pbccell(lengths...)
    end
    Rk = get(model, :Rk) do
        wrappos!(inv(cell) * model.R)
    end
    elements = get(model, :elements) do
        guesselements!(similar(model, String), model)
    end
    guesstopology!(topology, elements, Rk, cell; tol = tol)
end

function guesstopology!(topology::AbstractGraph,
        elements::AbstractVector{<:AbstractString}, Rk::AbstractMatrix{<:Real},
        cell::PBCCell; tol::Real = 0.1)
    dmax = (1.0+tol) * 2.0 * maximum(values(covalent_radius))
    g3 = append!(KspaceGrid3{Int}(cell, dmax), Rk, 1:ncols(Rk))
    Tk = similar(Rk, 3)
    T = similar(Rk, 3)
    J = Int[]
    for i = 1:ncols(Rk)
        Ak = @view Rk[:,i]
        for j in findnear!(J, g3, Ak)
            if i < j
                ri = covalent_radius[elements[i]]
                rj = covalent_radius[elements[j]]
                Bk = @view Rk[:,j]
                mintrans!(Tk, Ak, Bk)
                mul!(T, cell, Tk)
                if norm(T) < (1.0+tol) * (ri+rj)
                    pair!(topology, (i, j))
                end
            end
        end
    end
    topology
end
