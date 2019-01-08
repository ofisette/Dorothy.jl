# Properties of particles in molecular models

function wrapid(id::Integer, max::Integer)
	while id > max
		id -= max + 1
	end
	id
end

function unwrapids!(ids::AbstractVector{<:Integer}, max::Integer)
	offset = 0
	atzero = false
	for i in eachindex(ids)
		if ids[i] == 0
			if ! atzero
				offset += max + 1
			end
			atzero = true
		else
			atzero = false
		end
		ids[i] += offset
	end
	ids
end

function unwrapnames!(names::AbstractVector{<:AbstractString})
	for i in eachindex(names)
		name = names[i]
		if length(name) > 0 && isdigit(first(name))
			names[i] = chop(name, head=1, tail=0) * first(name)
		end
	end
	names
end

function namematcher(pattern::AbstractString)
	Istar = findall(==('*'), pattern)
	nstars = length(Istar)
	if nstars == 0
		==(pattern)
	elseif nstars == 1
		if length(pattern) == 1
			name -> true
		else
			if Istar[1] == 1
				suffix = chop(pattern, head=1, tail=0)
				name -> endswith(name, suffix)
			elseif Istar[1] == length(pattern)
				prefix = chop(pattern)
				name -> startswith(name, prefix)
			else
				error("wildcard * may only appear at start or end")
			end
		end
	elseif nstars == 2
		if Istar[1] != 1 || Istar[2] != length(pattern)
			error("wildcard * may only appear at start or end")
		else
			name -> occursin(name, pattern)
		end
	else
		error("wildcard * may only appear at start or end")
	end
end

function namematcher(s0::AbstractString, S::AbstractString...)
	F = [namematcher(s) for s in [s0, S...]]
	function (name)
		for f in F
			if f(name)
				return true
			end
		end
		false
	end
end

const hydrogen_name_pattern = "H*"

const vsite_name_pattern = ["MC*", "MN*", "MTRP*", "MW", "LP"]

const water_resname_pattern = ["HOH", "SOL", "TIP*", "WAT*"]

const acid_protein_resname_pattern = ["ASP", "GLU"]

const basic_protein_resname_pattern = ["ARG", "LYS"]

const charged_protein_resname_pattern =
		[acid_protein_resname_pattern..., basic_protein_resname_pattern...]

const polar_protein_resname_pattern =
		["ASN", "CYS", "GLN", "HIS", "SER", "THR", "TRP", "TYR"]

const hydrophobic_protein_resname_pattern =
		["ALA", "GLY", "ILE", "LEU", "MET", "PHE", "PRO", "VAL"]

const protein_resname_pattern =
		[charged_protein_resname_pattern...,
		 polar_protein_resname_pattern...,
		 hydrophobic_protein_resname_pattern...]

# J.L. King, T.H. Jukes. Science, 164(3881):788-98, 1969.
# doi:10.1126/science.164.3881.788
const vertebrate_aa_frequencies = Dict(
		"ALA" => 0.74,
		"ARG" => 0.42,
		"ASN" => 0.44,
		"ASP" => 0.59,
		"CYS" => 0.33,
		"GLU" => 0.58,
		"GLN" => 0.37,
		"GLY" => 0.74,
		"HIS" => 0.29,
		"ILE" => 0.38,
		"LEU" => 0.76,
		"LYS" => 0.72,
		"MET" => 0.18,
		"PHE" => 0.40,
		"PRO" => 0.50,
		"SER" => 0.81,
		"THR" => 0.62,
		"TRP" => 0.13,
		"TYR" => 0.33,
		"VAL" => 0.68)

const mainchain_name_pattern = ["N", "H", "CA", "C", "O", "OC1", "OC2", "OXT"]

const backbone_name_pattern = ["N", "CA", "C"]

const nuclacid_resname_pattern = ["A", "C", "G", "U", "DA", "DC", "DG", "DT"]

const lipid_resname_pattern = ["D3PC", "DLPC", "DLPE", "DLPG", "DMPA", "DMPC",
		"DMPG", "DOPC", "DOPE", "DOPG", "DOPS", "DPPC", "DPPE", "DPPG", "DPPS",
		"DSPG", "POPC", "POPE", "POPG", "POPS"]

const monatomic_ion_resname_pattern = ["NA", "CL", "CA", "MG", "K", "RB", "CS",
		"LI", "ZN", "H", "SR", "BA", "AL", "AG", "FE", "CU", "F", "BR", "I",
		"O", "S", "N", "P"]

const polyatomic_ion_resname_pattern =
		["CN", "CO3", "NO2", "NO3", "O2", "OH", "PO3", "PO4", "SO3", "SO4"]

const ion_resname_pattern =
		[monatomic_ion_resname_pattern..., polyatomic_ion_resname_pattern...]

const alphahelix_ss_pattern = "H"
const helix310_ss_pattern = "G"
const pihelix_ss_pattern = "I"
const turn_ss_pattern = "T"
const strand_ss_pattern = "E"
const bridge_ss_pattern = "B"
const coil_ss_pattern = "C"
const bend_ss_pattern = "S"

const helix_ss_pattern = [alphahelix_ss_pattern, helix310_ss_pattern,
		pihelix_ss_pattern, turn_ss_pattern]
const sheet_ss_pattern = [strand_ss_pattern, bridge_ss_pattern]
const loop_ss_pattern = [coil_ss_pattern, bend_ss_pattern]

const ishydrogen = namematcher(hydrogen_name_pattern)

isheavy(name) = !ishydrogen(name)

const isvsite = namematcher(vsite_name_pattern...)

const iswater = namematcher(water_resname_pattern...)

const isprotein = namematcher(sort(protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const isacidresidue = namematcher(sort(acid_protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const isbasicresidue = namematcher(sort(basic_protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const ischargedresidue = namematcher(sort(charged_protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const ispolarresidue = namematcher(sort(polar_protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const ishydrophobicresidue =
		namematcher(sort(hydrophobic_protein_resname_pattern,
		by=(resname -> vertebrate_aa_frequencies[resname]))...)

const ismainchainname = namematcher(mainchain_name_pattern...)

ismainchain(name, resname) = ismainchainname(name) && isprotein(resname)

issidechain(name, resname) = !ismainchainname(name) && isprotein(resname)

const isbackbonename = namematcher(backbone_name_pattern...)

isbackbone(name, resname) = isbackbonename(name) && isprotein(resname)

const isnuclacid = namematcher(nuclacid_resname_pattern...)

const islipid = namematcher(lipid_resname_pattern...)

const ision = namematcher(ion_resname_pattern...)

const ismonatomicion = namematcher(monatomic_ion_resname_pattern...)

const ispolyatomicion = namematcher(polyatomic_ion_resname_pattern...)

const ishelix = namematcher(helix_ss_pattern...)

const isalphahelix = namematcher(alphahelix_ss_pattern)

const ishelix310 = namematcher(helix310_ss_pattern)

const ispihelix = namematcher(pihelix_ss_pattern)

const isturn = namematcher(turn_ss_pattern)

const issheet = namematcher(sheet_ss_pattern...)

const isstrand = namematcher(strand_ss_pattern)

const isbridge = namematcher(bridge_ss_pattern)

const isloop = namematcher(loop_ss_pattern...)

const iscoil = namematcher(coil_ss_pattern)

const isbend = namematcher(bend_ss_pattern)

# J. Meija et al. Pure Appl. Chem., 88(3):265-91, 2016.
# doi:10.1515/pac-2015-0305
const standard_atomic_weights = Dict(
		"H"  => 1.008,
		"He" => 4.002,
		"Li" => 6.94,
		"Be" => 9.012,
		"B"  => 10.81,
		"C"  => 12.01,
		"N"  => 14.00,
		"O"  => 15.99,
		"F"  => 18.99,
		"Ne" => 20.18,
		"Na" => 22.99,
		"Mg" => 24.30,
		"Al" => 26.98,
		"Si" => 28.08,
		"P"  => 30.97,
		"S"  => 32.06,
		"Cl" => 35.45,
		"Ar" => 39.94,
		"K"  => 39.09,
		"Ca" => 40.07,
		"Sc" => 44.95,
		"Ti" => 47.86,
		"V"  => 50.94,
		"Cr" => 51.99,
		"Mn" => 54.93,
		"Fe" => 55.84,
		"Co" => 58.93,
		"Ni" => 58.69,
		"Cu" => 63.54,
		"Zn" => 65.38,
		"Ga" => 69.72,
		"Ge" => 72.63,
		"As" => 74.92,
		"Se" => 78.97,
		"Br" => 79.90,
		"Kr" => 83.79,
		"Rb" => 85.46,
		"Sr" => 87.62,
		"Y"  => 88.90,
		"Zr" => 91.22,
		"Nb" => 92.90,
		"Mo" => 95.95,
		"Ru" => 101.0,
		"Rh" => 102.9,
		"Pd" => 106.4,
		"Ag" => 107.8,
		"Cd" => 112.4,
		"In" => 114.8,
		"Sn" => 118.7,
		"Sb" => 121.7,
		"Te" => 127.6,
		"I"  => 126.9,
		"Xe" => 131.2,
		"Cs" => 132.9,
		"Ba" => 137.3,
		"La" => 138.9,
		"Ce" => 140.1,
		"Pr" => 140.9,
		"Nd" => 144.2,
		"Sm" => 150.3,
		"Eu" => 151.9,
		"Gd" => 157.2,
		"Tb" => 158.9,
		"Dy" => 162.5,
		"Ho" => 164.9,
		"Er" => 167.2,
		"Tm" => 168.9,
		"Yb" => 173.0,
		"Lu" => 174.9,
		"Hf" => 178.4,
		"Ta" => 180.9,
		"W"  => 183.8,
		"Re" => 186.2,
		"Os" => 190.2,
		"Ir" => 192.2,
		"Pt" => 195.0,
		"Au" => 196.9,
		"Hg" => 200.5,
		"Tl" => 204.3,
		"Pb" => 207.2,
		"Bi" => 208.9,
		"Th" => 232.0,
		"Pa" => 231.0,
		"U"  => 238.03)

# B. Cordero et al. Dalton Trans., (21):2832-8, 2008. doi:0.1039/b801115j
const covalent_radii = Dict(
	"H"  => 0.31,
	"He" => 0.28,
	"Li" => 1.28,
	"Be" => 0.96,
	"B"  => 0.84,
	"C"  => 0.76,
	"N"  => 0.71,
	"O"  => 0.66,
	"F"  => 0.57,
	"Ne" => 0.58,
	"Na" => 1.66,
	"Mg" => 1.41,
	"Al" => 1.21,
	"Si" => 1.11,
	"P"  => 1.07,
	"S"  => 1.05,
	"Cl" => 1.02,
	"Ar" => 1.06,
	"K"  => 2.03,
	"Ca" => 1.76,
	"Sc" => 1.70,
	"Ti" => 1.60,
	"V"  => 1.53,
	"Cr" => 1.39,
	"Mn" => 1.39,
	"Fe" => 1.32,
	"Co" => 1.26,
	"Ni" => 1.24,
	"Cu" => 1.32,
	"Zn" => 1.22,
	"Ga" => 1.22,
	"Ge" => 1.20,
	"As" => 1.19,
	"Se" => 1.20,
	"Br" => 1.20,
	"Kr" => 1.16,
	"Rb" => 2.20,
	"Sr" => 1.95,
	"Y"  => 1.90,
	"Zr" => 1.75,
	"Nb" => 1.64,
	"Mo" => 1.54,
	"Tc" => 1.47,
	"Ru" => 1.46,
	"Rh" => 1.42,
	"Pd" => 1.39,
	"Ag" => 1.45,
	"Cd" => 1.44,
	"In" => 1.42,
	"Sn" => 1.39,
	"Sb" => 1.39,
	"Te" => 1.38,
	"I"  => 1.39,
	"Xe" => 1.40,
	"Cs" => 2.44,
	"Ba" => 2.15,
	"La" => 2.07,
	"Ce" => 2.04,
	"Pr" => 2.03,
	"Nd" => 2.01,
	"Pm" => 1.99,
	"Sm" => 1.98,
	"Eu" => 1.98,
	"Gd" => 1.96,
	"Tb" => 1.94,
	"Dy" => 1.92,
	"Ho" => 1.92,
	"Er" => 1.89,
	"Tm" => 1.90,
	"Yb" => 1.87,
	"Lu" => 1.87,
	"Hf" => 1.75,
	"Ta" => 1.70,
	"W"  => 1.62,
	"Re" => 1.51,
	"Os" => 1.44,
	"Ir" => 1.41,
	"Pt" => 1.36,
	"Au" => 1.36,
	"Hg" => 1.32,
	"Tl" => 1.45,
	"Pb" => 1.46,
	"Bi" => 1.48,
	"Po" => 1.40,
	"At" => 1.50,
	"Rn" => 1.50,
	"Fr" => 2.60,
	"Ra" => 2.21,
	"Ac" => 2.15,
	"Th" => 2.06,
	"Pa" => 2.00,
	"U"  => 1.96,
	"Np" => 1.90,
	"Pu" => 1.87,
	"Am" => 1.80,
	"Cm" => 1.69)
