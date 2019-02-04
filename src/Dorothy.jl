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
		SelectionCache, SelectionMode, Selector, Selectors, DorothyIO

include("utils.jl")
include("properties.jl")
include("headers.jl")
include("graphs.jl")
include("geometry.jl")
include("pbc.jl")
include("multicolls.jl")

using .Utils
using .Properties
using .Headers
using .Graphs
using .Geometry
using .PBC
using .Multicollections

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

function __init__()
	register_gromos87()
	register_pdb()
	register_trr()
	register_xtc()
	register_ndx()
end

end # module
