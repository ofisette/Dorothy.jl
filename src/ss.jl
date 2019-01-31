module SS

using ..Dorothy
using ..Dorothy.Geometry
using ..Dorothy.PDB
using ..Dorothy.Selectors
using ..Dorothy.Utils

export guessss!

abstract type SSGuessStrategy end

guessss!(model::ParticleCollection) = guessss!(get!(model, :SS, ""), model)

guessss!(model::ParticleCollection, strategy::SSGuessStrategy) =
		guessss!(get!(model, :SS, ""), model, strategy)

struct SSByStride <: SSGuessStrategy
	stridepath::String

	SSByStride(stridepath::AbstractString = "stride") = new(stridepath)
end

guessss!(SS::AbstractVector{<:AbstractString}, model::ParticleCollection) =
		guessss!(SS, model, SSByStride())

function guessss!(SS::AbstractVector{<:AbstractString},
		model::ParticleCollection, strategy::SSByStride)
	@boundscheck length(SS) == length(model) ||
			error("size mismatch between model and output array")
	I = findall(map(Protein, model))
	protein = view(model, I)
	chainids = get(protein, :chainids, Repeat("", length(protein)))
	residues = mcrp(protein).flath2.tree
	stride!(view(SS, I), protein.names, protein.resids, protein.resnames,
			chainids, protein.R, residues, strategy.stridepath)
end

function stride!(SS::AbstractVector{<:AbstractString},
		names::AbstractVector{<:AbstractString},
		resids::AbstractVector{<:Integer},
		resnames::AbstractVector{<:AbstractString},
		chainids::AbstractVector{<:AbstractString},
		R::AbstractVector{Vector3D},
		residues::AbstractVector{<:AbstractVector{<:Integer}},
		stridepath::AbstractString)
	n = length(SS)
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
					ss = prevss
				else
					ss = tokens[6]
				end
				push!(striderecords, (resid, resname, ss))
				prevss = ss
			end
		end
	end
	nstride = length(striderecords)
	nres = length(residues)
	nstride == nres ||
			error("cannot assign $(nstride) stride records to $(nres) residues")
	for (I, (resid, resname, ss)) in zip(residues, striderecords)
		resids[I[1]] == resid ||
				error("stride resid $(resid) does not match $(resids[I[1]])")
		resnames[I[1]] == resname ||
				error("stride resname mismatch at resid $(resid[I[1]])")
		for i in I
			SS[i] = ss
		end
	end
	SS
end

end
