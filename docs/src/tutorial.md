# Tutorial

This page will give you a quick overview of the basic usage of Dorothy before
you read through the *Manual* or selected pages thereof. You will learn to read
molecular models from files, query and modify particle properties, select
subsets of the model, write a subset to a new file, and compute simple
geometrical properties from particle positions.

## Reading molecular models

We assume Dorothy was already added to your Julia environment (see
[Installation](@ref)). We load the necessary modules, then read a molecular
model, an MHC-I protein in solvent, from a PBD file distributed with the
package.

```julia-repl
julia> using Dorothy, Formats

julia> model = readf("$(Dorothy.datapath)/MHC.pdb")
68336-particle 9-key MolecularModel:
 bfactors
 chainids
 elements
 ids
 names
 occupancies
 R
 resids
 resnames
```

There is quite a bit going on in these few lines. The `Dorothy` top-level module
provides basic data types, such as `MolecularModel`. The `Formats` module
provides `readf`, a generic function to detect and read formatted data. Here,
`readf` delegates reading `MHC.pdb` to the `Dorothy.PDB` submodule. The
resulting `MolecularModel` contains 68336 particles (atoms) with 9 properties
typical of PDB structures.

## Particle properties

Particle properties are arrays:

```julia-repl
julia> model.names
68336-element Array{String,1}:
 "MN1"
 "MN2"
 "N"  
 "H1"
 "H2"
 "H3"
 "CA"
 "HA1"
 "HA2"
 "C"  
 "O"  
 ⋮    
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"
 "CL"

julia> model.resnames[1]
"GLY"

julia> model.elements[1:10]
10-element Array{String,1}:
 ""
 ""
 "N"
 "H"
 "H"
 "H"
 "C"
 "H"
 "H"
 "C"
```

Individual particles can also be queried by indexing the model. Properties
become scalar quantities and use the singular form rather the plural.

```julia-repl
julia> p7 = model[7]
9-key MolecularModel particle:
 bfactor
 chainid
 element
 id
 name
 occupancy
 R
 resid
 resname

julia> p7.name
"CA"
```

Properties can be set on the model or individual particles.

```julia-repl
julia> p7.occupancy = 0.0
0.0

julia> model.resids[1:11] .= 0
11-element view(::Array{Int64,1}, 1:11) with eltype Int64:
 0
 0
 0
 0
 0
 0
 0
 0
 0
 0
 0
```

Particle positions are represented by `Vector3D` objects from the `Geometry`
submodule. They are fixed-size, immutable, convertible to other `AbstractVector`
types, and have `x`, `y` and `z` fields in addition to array-like indexing.

```julia-repl
julia> using
julia> model.R
68336-element Array{Dorothy.Geometry.Vector3D,1}:
 [88.98, 81.17, 48.96]
 [88.5, 81.43, 48.36]  
 [88.79, 81.31, 48.62]
 [89.09, 80.96, 49.53]
 [88.32, 82.19, 48.74]
 [88.16, 80.64, 48.2]  
 [89.96, 81.48, 47.76]
 [90.72, 82.16, 48.24]
 [90.48, 80.5, 47.58]  
 [89.49, 82.08, 46.44]
 [89.17, 83.26, 46.35]
 ⋮                     
 [46.27, 45.06, 37.91]
 [48.76, 62.74, 3.36]  
 [58.17, 32.9, 6.03]   
 [82.55, 11.78, 63.57]
 [57.78, 8.74, 39.17]  
 [32.68, 17.64, 6.86]  
 [51.55, 60.28, 0.5]   
 [31.73, 14.62, 24.97]
 [113.71, 82.34, 42.32]
 [82.43, 65.82, 68.69]
 [110.61, 81.62, 38.56]

julia> model[1].R
3-element Dorothy.Geometry.Vector3D:
 88.98
 81.17
 48.96

julia> model[1].R[1]
88.98

julia> model[1].R.x
88.98

julia> model[1].R = Vector3D(1.0, 1.0, 1.0)
3-element Vector3D:
 1.0
 1.0
 1.0

julia> model[1].R = [1.0, 1.0, 1.0]
3-element Array{Float64,1}:
 1.0
 1.0
 1.0

julia> model[1].R.x = 3.0
ERROR: setfield! immutable struct of type Vector3D cannot be changed
[...]
```

## Model headers

In addition to particle properties, each model has a header that stores global
properties, such as the unit cell (for periodic systems), a title, etc. This
`MolecularModelHeader` behaves like a dictionary, but also allows named property
access and enforces appropriate types for properties.

```julia-repl
julia> model.header
2-entry MolecularModelHeader:
 cell
 title

julia> model.header.title
"MHC-I protein, water and ions"

julia> model.header.title = "Immunity-related protein in solvent box"
"Immunity-related protein in solvent box"

julia> model.header.title = 42
ERROR: MethodError: no method matching headerval(::MolecularModelHeader, ::Val{:title}, ::Int64)
[...]
```

## Particle selections

Selecting a subset of a molecular model to perform calculations on this subset
is a very common operation. The `Base.view` function, which returns a view into a parent array, can also be used with particle collections.

```julia-repl
julia> prot = view(model, 1:6539)
6539-particle 9-key MolecularModel view:
 bfactors
 chainids
 elements
 ids
 names
 occupancies
 R
 resids
 resnames

julia> protC = view(prot, 6391:6539)
149-particle 9-key MolecularModel view:
 bfactors
 chainids
 elements
 ids
 names
 occupancies
 R
 resids
 resnames
```

Here, we create two views, the first for all protein atoms in the model, the
second only for the third polypeptide chain in the model. Typically, however,
you do not know in advance the particle indices to target. Instead, you wish to
perform selections such as “all protein atoms”, or “all atoms from chain C”. To
do this, the `view` function as it applies to particle collections is extended
to accept `Selector` objects which query particle properties to determine if
they should be included in a selection. The `Selectors` submodule contains
predefined selectors.

```julia-repl
julia> using Dorothy.Selectors

julia> prot == view(model, Protein)
true

julia> protC == view(prot, ChainId("C"))
true
```

As you can see, the above two views can be created using selections in a much
more convenient manner. Selectors can be combined to perform powerful operations
on molecular models.

```julia-repl
julia> selector = Water & Restrict(Within(5.0, of=Protein), by=Residue)
[...]

julia> length(view(model, selector)) / 3
1921.0
```

There are 1921 water molecules fully enclosed within 5.0 Å of any protein atom.

## Writing molecular models

Particle collections (both models and views) can be written to files using
`writef` from the `Formats` module. Make sure to use a writable location for the
output file in the next code example.

```julia-repl
julia> writef("MHC_protC.pdb", protC)
12183
```

## Geometry calculations
