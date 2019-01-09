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
		ParticleCollection, Particle, MolecularTrajectory,

		hierarchy, eachchain, eachresidue, eachlocalresidue, eachfragment,
		chainat, residueat, fragmentat,

		wrapid, unwrapids!, unwrapnames!,

		namematcher, ishydrogen, isheavy, isvsite, iswater, isprotein,
		isacidresidue, isbasicresidue, ischargedresidue, ispolarresidue,
		ishydrophobicresidue, ismainchain, issidechain, isbackbone, isnuclacid,
		islipid, ision, ismonatomicion, ispolyatomicion, isalphahelix,
		ishelix310, ispihelix, isturn, isstrand, isbridge, iscoil, isbend,
		ishelix, issheet, isloop,

		guesselements!, guessmasses!, guessss!, guesstopology!,

		SelectionCache, emptyframe!, SelectBy, Selector, pick, Selectors,

		DorothyIO

include("utils.jl")
include("header.jl")
include("graphs.jl")
include("geometry.jl")
include("pbc.jl")
include("multicolls.jl")

using .Utils
using .Headers
using .Graphs
using .Geometry
using .PBC
using .Multicollections

include("models.jl")
include("properties.jl")
include("guessing.jl")
include("selection.jl")

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

function __init__()
	register_gromos87()
	register_pdb()
	register_trr()
	register_xtc()
	register_ndx()
end

end # module
