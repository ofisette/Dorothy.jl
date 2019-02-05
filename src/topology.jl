module Topology

using ..Dorothy
using ..Dorothy.Atoms
using ..Dorothy.Geometry
using ..Dorothy.Graphs
using ..Dorothy.PBC
using ..Dorothy.Properties
using ..Dorothy.Selectors
using ..Dorothy.Utils

export
		TopologyInferenceStrategy, infertopology, infertopology!,
		TopologyFromCovalentRadii

abstract type TopologyInferenceStrategy end

infertopology(model::ParticleCollection) =
		infertopology!(Graph(length(model)), model)

infertopology(model::ParticleCollection, strategy::TopologyInferenceStrategy) =
		infertopology!(Graph(length(model)), model, strategy)

infertopology!(model::ParticleCollection) =
		infertopology!(get!(model, :topology, []), model)

infertopology!(model::ParticleCollection, strategy::TopologyInferenceStrategy) =
		infertopology!(get!(model, :topology, []), model, strategy)

struct TopologyFromCovalentRadii <: TopologyInferenceStrategy
	radii::Dict{String,Float64}
	rdef::Float64
	rtol::Float64
	dmaxXX::Float64
	dmaxXH::Float64
	dmaxHH::Float64

	function TopologyFromCovalentRadii(;
			radii::AbstractDict{<:AbstractString,<:Real} = covalent_radii,
			rdef::Real = radii["C"], tol::Real = 0.2)
		@boundscheck tol >= 0.0 || error("expected positive tolerance")
		rtol = 1.0 + tol
		rmax = maximum(values(radii))
		rH = radii["H"]
		dmaxXX = rtol * (rmax + rmax)
		dmaxXH = rtol * (rmax + rH)
		dmaxHH = rtol * (rH + rH)
		new(radii, rdef, rtol, dmaxXX, dmaxXH, dmaxHH)
	end
end

infertopology!(topology::AbstractGraph, model::ParticleCollection) =
		infertopology!(topology, model, TopologyFromCovalentRadii())

function infertopology!(topology::AbstractGraph, model::ParticleCollection,
		strategy::TopologyFromCovalentRadii)
	n = length(topology)
	@boundscheck length(model) == n ||
			error("size mismatch between model and output array")
	elements = infermissingelements(model)
	radii = [get(strategy.radii, elements[i], strategy.rdef) for i = 1:n]
	H = map(Hydrogen, model)
	X = .!H
	cell, Rw, Kw = pbcstrategy(model)
	XXlattice = proxilattice(Kw, cell, strategy.dmaxXX)
	XHlattice = proxilattice(Kw, cell, strategy.dmaxXH)
	HHlattice = proxilattice(Kw, cell, strategy.dmaxHH)
	topcov!(topology, radii, Rw, Kw, cell, XXlattice, X, X, strategy.rtol)
	topcov!(topology, radii, Rw, Kw, cell, XHlattice, X, H, strategy.rtol)
	topcov!(topology, radii, Rw, Kw, cell, HHlattice, H, H, strategy.rtol)
	topology
end

function topcov!(topology::AbstractGraph, radii::AbstractVector{<:Real},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, lattice::ProximityLattice,
		Itargets::AbstractVector{Bool}, Jtargets::AbstractVector{Bool},
		rtol::Real)
	n = length(topology)
	@boundscheck begin
		length(radii) == length(Rw) == length(Kw) == length(lattice) ==
				length(Itargets) == length(Jtargets) == n ||
				error("size mismatch between property arrays")
		rtol >= 1.0 || error("expected >= 1.0 relative tolerance")
	end
	J = Int[]
	for i = 1:n
		if Itargets[i]
			for j in findnear!(J, lattice, i)
				if Jtargets[j] && (i != j)
					if mindist(Rw[i], Kw[i], Rw[j], Kw[j], cell) <
							rtol * (radii[i] + radii[j])
						pair!(topology, (i, j))
					end
				end
			end
		end
	end
end

end
