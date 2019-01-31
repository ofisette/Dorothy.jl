module Topology

using ..Dorothy
using ..Dorothy.Geometry
using ..Dorothy.Graphs
using ..Dorothy.PBC

export guesstopology!

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

end
