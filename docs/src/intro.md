# Introduction

Dorothy provides data structures and algorithms for the scripted or interactive
manipulation and analysis of molecular systems, with a particular emphasis on
biological macromolecules. You can use Dorothy to build or modify systems for
starting molecular dynamics simulations, to analyse trajectories output by such
simulations, or to examine structures from the Protein Data Bank.

## Main package features

- Array-like data structures representing particle collections
  - Access particle properties as arrays
  - Create views over particle subsets
  - Splice particles to add/remove them from a model
- Read and write molecular structures in common file formats
  - Protein Data Bank (`.pdb`)
  - Gromos87 (`.gro`)
- Read and write molecular trajectories in common file formats
  - Gromacs uncompressed trajectories (`.trr`)
  - Gromacs compressed trajectories (`.xtc`, read-only)
- Automatic format detection with support for compressed files
- Hierarchical iterators over models
  - Model/chain/residue/particle
  - Model/fragment/particle
- Particle selection
  - Property-based selectors (name, id, mass, secondary structure, etc.)
  - Hierarchy-based selectors (complete or partial chain, residue, etc.)
  - Distance-based selectors
  - Logical selectors
  - Predefined convenience selectors (protein, N/C-termini, solvent, etc.)
- Geometry
  - Distances, angles, dihedrals, centroids
  - Transformations: translation, rotation, scaling
  - RMSD and superposition
  - Line and plane fitting
- Triclinic periodic boundary conditions
  - Minimum particle-particle distance across periodic images
  - Wrap and unwrap particles into/from the unit cell
- Efficient distance calculations using neighbor lists and cell grids
- Automated topology assignment
- Automated secondary structure assignment (using STRIDE)

## Using this documentation

The documentation for Dorothy is split in three main parts. The first part,
*Manual*, describes all the important features of the package, by topic and in a
logical order. These pages form the core of the documentation and are especially
recommended for new users. Before diving into the manual, however, it is
suggested to go through the short [Tutorial](@ref) that follows this
introduction.

The second part, *Development*, discusses the general philosophy behind Dorothy
and how to extend the package by implementing new features and writing generic,
reusable code. Most users can ignore this.

The third part, *Reference*, covers the public interface of Dorothy and its
submodules. The same documentation is available from the REPL using Julia’s help
mechanism; `?Dorothy` can be used as a starting point to navigate through the
submodules and their public variables.
