"""
Molecular structure manipulation and analysis
"""
module Dorothy

using LinearAlgebra
using Statistics
using Formats
using FormatStreams

export
		MolecularModelHeader, MolecularModel, MolecularModelView,
		ParticleCollection, Particle, MolecularTrajectory, pbcstrategy,
		pbcstrategy!, SelectionCache, SelectionMode, Selector, Selectors,
		DorothyIO, @DorothyAll

include("utils.jl")
include("properties.jl")
include("headers.jl")
include("graphs.jl")
include("multicolls.jl")
include("geometry.jl")
include("pbc.jl")
include("neighbors.jl")

using .Utils
using .Properties
using .Headers
using .Graphs
using .Multicollections
using .Geometry
using .PBC
using .Neighbors

include("models.jl")
include("atoms.jl")
include("hierarchies.jl")

using .Atoms
using .Hierarchies

include("selections.jl")

using .Selectors

struct DorothyIO <: Formats.FormatHandler end

include("formats/gromos87.jl")
include("formats/pdb.jl")
include("formats/xdr.jl")
include("formats/trr.jl")
include("formats/xtc.jl")
include("formats/ndx.jl")

using .Gromos87
using .PDB
using .TRR
using .XTC
using .NDX

include("ss.jl")
include("topology.jl")

using .SecondaryStructure
using .Topology

macro DorothyAll()
	return :( using LinearAlgebra,
			Statistics,
			Formats,
			FormatCodecs,
			FormatStreams,
			Dorothy,
			Dorothy.Utils,
			Dorothy.Properties,
			Dorothy.Graphs,
			Dorothy.Multicollections,
			Dorothy.Geometry,
			Dorothy.PBC,
			Dorothy.Neighbors,
			Dorothy.Atoms,
			Dorothy.Hierarchies,
			Dorothy.Selectors,
			Dorothy.Gromos87,
			Dorothy.PDB,
			Dorothy.TRR,
			Dorothy.XTC,
			Dorothy.NDX,
			Dorothy.SecondaryStructure,
			Dorothy.Topology )
end

datapath = joinpath(@__DIR__, "..", "data")

function __init__()
	register_gromos87()
	register_pdb()
	register_trr()
	register_xtc()
	register_ndx()
end

end # module
