# Protein secondary structure

# SS symbols
#   Helix
#     H: alpha helix
#     G: 3_10 helix
#     I: pi helix
#     T: turn
#   Sheet
#     E: strand
#     B: bridge
#   Loop
#     C: coil
#     S: bend

module SecondaryStructure

using ..Dorothy
using ..Dorothy.Geometry
using ..Dorothy.Hierarchies
using ..Dorothy.PDB
using ..Dorothy.Selectors
using ..Dorothy.Utils

export SSInferenceStrategy, inferss, inferss!, SSByStride

abstract type SSInferenceStrategy end

inferss(model::ParticleCollection) =
		inferss!(Vector{String}(undef, length(model)), model)

inferss!(model::ParticleCollection) = inferss!(get!(model, :SS, ""), model)

inferss!(model::ParticleCollection, strategy::SSInferenceStrategy) =
		inferss!(get!(model, :SS, ""), model, strategy)

struct SSByStride <: SSInferenceStrategy
	stridepath::String

	SSByStride(; stridepath::AbstractString = "stride") = new(stridepath)
end

inferss!(SS::AbstractVector{<:AbstractString}, model::ParticleCollection) =
		inferss!(SS, model, SSByStride())

function inferss!(SS::AbstractVector{<:AbstractString},
		model::ParticleCollection, strategy::SSByStride)
	@boundscheck length(SS) == length(model) ||
			error("size mismatch between model and output array")
	I = findall(map(Protein, model))
	protein = view(model, I)
	chainids = get(protein, :chainids, Repeated("", length(protein)))
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
				R, Repeated(1.0, n), Repeated(0.0, n), Repeated("", n),
				nothing, true, Ref(0))
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
