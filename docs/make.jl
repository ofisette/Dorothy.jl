using Documenter
using Dorothy

docargs = (
	sitename="Dorothy.jl",
	pages = [
		"Home" => "index.md",
		"intro.md",
		"tutorial.md",
		"Manual" => [
			"module.md",
			"models.md"
		],
		"Development" => [],
		"Reference" => [],
		"Appendices" => [
			"refs.md",
			"license.md"
		]
	]
)

makedocs(;
	docargs...,
	build="build/web"
)

makedocs(;
	docargs...,
	build="build/local",
	format = Documenter.HTML(prettyurls = false)
)
