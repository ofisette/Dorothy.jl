# Tutorial

This page will give you a quick overview of the basic usage of Dorothy before
you read through the *Manual* or selected pages thereof. You will learn to read
molecular models from files, query and modify particle properties, select
subsets of the model, write a subset to a new file, and compute simple
geometrical properties from particle positions.

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
