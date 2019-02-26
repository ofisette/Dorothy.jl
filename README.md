# Dorothy.jl

Julia package for molecular structure manipulation and analysis

## License

You can use Dorothy under the terms of the MIT License; see
[`LICENSE.md`](https://github.com/ofisette/Dorothy.jl/blob/master/LICENSE.md) in
the project files, or [License](https://ofisette.github.io/Dorothy.jl/license/)
in the online documentation.

## Status

Dorothy is in early development, and targets Julia 1.x.

## Installation

Dorothy is not a registered package. You can add it to your Julia environment by
giving the URL to its repository. The unregistered packages Formats,
FormatCodecs and FormatStreams are required as dependencies and can be installed
in the same manner.

```julia
using Pkg

Pkg.add("https://github.com/ofisette/Formats.jl")
Pkg.add("https://github.com/ofisette/FormatCodecs.jl")
Pkg.add("https://github.com/ofisette/FormatStreams.jl")
Pkg.add("https://github.com/ofisette/Dorothy.jl")
```

## Documentation

Documentation is available [online](https://ofisette.github.io/Dorothy.jl/), and
can be built from the project files. Documentation for the the public API is
also available directly from the Julia REPL.

## Citing

If you use Dorothy in your research, please cite the software in your
publications. There is currently no published work describing Dorothy, but it
can still be cited as follows:

Olivier Fisette. *Dorothy: Julia package for molecular structure manipulation
and analysis*, version 0.1, 2019. `https://github.com/ofisette/Dorothy.jl`

## Community

Dorothy is developed by [Olivier Fisette](mailto:olivier.fisette@rub.de) in the
[Molecular Simulation Group](https://molecular-simulation.org/) of Lars V.
Schäfer at the Center for Theoretical Chemistry of Ruhr-University Bochum,
Germany.

Contributions, bug reports and feature suggestions are welcome. Development is
tracked in the project’s [repository](https://github.com/ofisette/Dorothy.jl).
