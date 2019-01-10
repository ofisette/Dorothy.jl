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

function guesselements!(elements::AbstractVector{<:AbstractString},
		model::ParticleCollection)
	@boundscheck length(elements) == length(model) ||
			error("size mismatch between model and output array")
	guesselements!(elements, model.names, model.resnames)
end

function guesselements!(elements::AbstractVector{<:AbstractString},
		names::AbstractVector{<:AbstractString},
		resnames::AbstractVector{<:AbstractString})
	@boundscheck begin
		length(elements) == length(names) == length(resnames) ||
				error("size mismatch between property arrays")
	end
	for i in eachindex(elements)
		elements[i] = guesselement(names[i], resnames[i])
	end
	elements
end

function guessmass(name::AbstractString, element::AbstractString)
	get(standard_atomic_weights, element) do
		if isvsite(name)
			0.0
		else
			error("could not guess mass from name or element")
		end
	end
end

guessmasses!(model::ParticleCollection) =
		guessmasses!(get!(model, :masses, undef), model)

function guessmasses!(masses::AbstractVector{<:Real}, model::ParticleCollection)
	@boundscheck length(masses) == length(model) ||
			error("size mismatch between model and output array")
	guessmasses!(masses, model.names,
			get(model, :elements, guesselements!(model)))
end

function guessmasses!(masses::AbstractVector{<:Real},
		names::AbstractVector{<:AbstractString},
		elements::AbstractVector{<:AbstractString})
	@boundscheck begin
		length(masses) == length(names) == length(elements) ||
				error("size mismatch between property arrays")
	end
	for i in eachindex(masses)
		masses[i] = guessmass(names[i], elements[i])
	end
	masses
end

abstract type SSGuessStrategy end

guessss!(model::ParticleCollection) = guessss!(get!(model, :SS, ""), model)

guessss!(model::ParticleCollection, strategy::SSGuessStrategy) =
		guessss!(get!(model, :SS, ""), model, strategy)

struct SSByStride <: SSGuessStrategy
	stridepath::String

	SSByStride(stridepath::AbstractString = "stride") = new(stridepath)
end

guessss!(ss::AbstractVector{<:AbstractString}, model::ParticleCollection) =
		guessss!(ss, model, SSByStride())

function guessss!(ss::AbstractVector{<:AbstractString},
		model::ParticleCollection, strategy::SSByStride)
	@boundscheck length(ss) == length(model) ||
			error("size mismatch between model and output array")
	I = findall(map(Protein, model))
	protein = view(model, I)
	chainids = get(protein, :chainids, Repeat("", length(protein)))
	Ires2ps = eachresidue(protein).Igroup2items
	stride!(view(ss, I), protein.names, protein.resids, protein.resnames,
			chainids, protein.R, Ires2ps, strategy.stridepath)
end

function stride!(ss::AbstractVector{<:AbstractString},
		names::AbstractVector{<:AbstractString},
		resids::AbstractVector{<:Integer},
		resnames::AbstractVector{<:AbstractString},
		chainids::AbstractVector{<:AbstractString},
		R::AbstractVector{Vector3D},
		Ires2ps::AbstractVector{<:AbstractVector{<:Integer}},
		stridepath::AbstractString)
	n = length(ss)
	@boundscheck length(names) == length(resids) == length(resnames) ==
			length(chainids) == length(R) == n ||
			error("size mismatch between property arrays")
	striderecords = Tuple{Int,String,String}[]
	mktemp() do pdbpath, io
		writepdb(io, nothing, nothing, 1:n, names, resids, resnames, chainids,
				R, Repeat(1.0, n), Repeat(0.0, n), Repeat("", n), nothing, true,
				Ref(0))
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
	nres = length(Ires2ps)
	nstride == nres ||
			error("cannot assign $(nstride) stride records to $(nres) residues")
	for (I, (resid, resname, ssi)) in zip(Ires2ps, striderecords)
		resids[I[1]] == resid ||
				error("stride resid $(resid) does not match $(resids[I[1]])")
		resnames[I[1]] == resname ||
				error("stride resname mismatch at resid $(resid[I[1]])")
		for i in I
			ss[i] = ssi
		end
	end
	ss
end

abstract type TopologyGuessStrategy end

guesstopology!(model::ParticleCollection) =
		guesstopology!(get!(model, :topology, []), model)

guesstopology!(model::ParticleCollection, strategy::TopologyGuessStrategy) =
		guesstopology!(get!(model, :topology, []), model, strategy)

struct TopologyByCovalentRadius
	radii::Dict{String,Float64}
	tol::Float64
	dmax::Float64

	function TopologyByCovalentRadius(
			radii::AbstractDict{<:AbstractString,<:Real} = covalent_radii,
			tol::Real = 0.1)
		@boundscheck tol > 0.0 || error("expected strictly positive tolerance")
		dmax = (1.0+tol) * 2.0 * maximum(values(radii))
		new(radii, tol, dmax)
	end
end

guesstopology!(topology::AbstractGraph, model::ParticleCollection) =
		guesstopology!(topology, model, TopologyByCovalentRadius())

function guesstopology!(topology::AbstractGraph, model::ParticleCollection,
		strategy::TopologyByCovalentRadius)
	@boundscheck length(topology) == length(model) ||
			error("size mismatch between model and output array")
	elements = get(model, :elements) do
		guesselements!(similar(model, String), model)
	end
	cell = get(model.header, :cell, nothing)
	Rw, Kw = pbcpos(model.R, cell)
	pg = posgrid(Rw, Kw, cell, strategy.dmax)
	topcov!(topology, elements, Rw, Kw, cell, pg, strategy.tol, strategy.radii)
end

function topcov!(topology::AbstractGraph,
		elements::AbstractVector{<:AbstractString},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, pg::PositionGrid, tol::Real,
		radii::AbstractDict{<:AbstractString,<:Real})
	@boundscheck length(topology) == length(elements) == length(Rw) ==
			length(Kw) || error("size mismatch between property arrays")
	J = Int[]
	for i in eachindex(Rw)
		for j in findnear!(J, pg, i)
			if i != j && elements[i] != "" && elements[j] != ""
				ri = radii[elements[i]]
				rj = radii[elements[j]]
				if mindist(Rw[i], Kw[i], Rw[j], Kw[j], cell) <
						(1.0+tol) * (ri+rj)
					pair!(topology, (i,j))
				end
			end
		end
	end
	topology
end
