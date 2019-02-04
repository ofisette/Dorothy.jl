module Topology

using ..Dorothy
using ..Dorothy.Atoms
using ..Dorothy.Geometry
using ..Dorothy.Graphs
using ..Dorothy.PBC
using ..Dorothy.Properties
using ..Dorothy.Utils

export TopologyInferenceStrategy, infertopology!, TopologyFromCovalentRadii

abstract type TopologyInferenceStrategy end

infertopology!(model::ParticleCollection) =
		infertopology!(get!(model, :topology, []), model)

infertopology!(model::ParticleCollection, strategy::TopologyInferenceStrategy) =
		infertopology!(get!(model, :topology, []), model, strategy)

struct TopologyFromCovalentRadii <: TopologyInferenceStrategy
	radii::Dict{String,Float64}
	rdef::Float64
	tol::Float64
	dmax::Float64

	function TopologyFromCovalentRadii(;
			radii::AbstractDict{<:AbstractString,<:Real} = covalent_radii,
			rdef::Real = radii["C"], tol::Real = 0.2)
		@boundscheck tol >= 0.0 || error("expected positive tolerance")
		dmax = (1.0+tol) * 2.0 * maximum(values(radii))
		new(radii, rdef, tol, dmax)
	end
end

infertopology!(topology::AbstractGraph, model::ParticleCollection) =
		infertopology!(topology, model, TopologyFromCovalentRadii())

function infertopology!(topology::AbstractGraph, model::ParticleCollection,
		strategy::TopologyFromCovalentRadii)
	@boundscheck length(topology) == length(model) ||
			error("size mismatch between model and output array")
	elements = infermissingelements(model)
	cell = pbccell(get(model.header, :cell, nothing))
	Rw, Kw = pbcpos(model.R, cell)
	lattice = proxilattice(Kw, cell, strategy.dmax)
	topcov!(topology, elements, Rw, Kw, cell, lattice, strategy.radii,
			strategy.rdef, strategy.tol)
end

function topcov!(topology::AbstractGraph,
		elements::AbstractVector{<:AbstractString},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, lattice::ProximityLattice,
		radii::AbstractDict{<:AbstractString,<:Real}, rdef::Real, tol::Real)
	@boundscheck length(topology) == length(elements) == length(Rw) ==
			length(Kw) || error("size mismatch between property arrays")
	tolf = 1.0 + tol
	J = Int[]
	for i in eachindex(Rw)
		for j in findnear!(J, lattice, i)
			if i != j && !((i, j) in topology)
				ri = get(radii, elements[i], rdef)
				rj = get(radii, elements[j], rdef)
				if mindist(Rw[i], Kw[i], Rw[j], Kw[j], cell) < tolf * (ri+rj)
					pair!(topology, (i, j))
				end
			end
		end
	end
	topology
end

end
