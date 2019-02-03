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

		mcrptree, mcrp, chains, residues, mfptree, mfp, fragments,
		chainat, residueat, fragmentat,

		standard_atomic_weights, covalent_radii, vertebrate_aa_frequencies,

		wrapid, unwrapids!, unwrapnames!,

		namematcher, isheavy, isvsite, iswater, isprotein, isacidresidue,
		isbasicresidue, ischargedresidue, ispolarresidue, ishydrophobicresidue,
		ismainchain, issidechain, isbackbone, isnuclacid, islipid, ision,
		ismonatomicion, ispolyatomicion, isalphahelix, ishelix310, ispihelix,
		isturn, isstrand, isbridge, iscoil, isbend, ishelix, issheet, isloop,

		inferelement, inferelements!, infermissingelements,
		infermissingelements!, infermass, infermasses!,

		SelectionCache, emptyframe!, SelectionMode, Selector, Selectors,

		DorothyIO

include("utils.jl")
include("headers.jl")
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
include("selections.jl")

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

using .SS
using .Topology

function __init__()
	register_gromos87()
	register_pdb()
	register_trr()
	register_xtc()
	register_ndx()
end

end # module
