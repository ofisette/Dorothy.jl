using Documenter
using Dorothy

docargs = (
	sitename="Dorothy.jl",
	pages = [
		"Home" => "index.md",
		"intro.md",
		"tutorial.md",
		"Manual" => [
			"models.md"
		],
		"Development" => [],
		"Appendices" => [
			"refs.md",
			"license.md"
		]
	],
	assets = ["assets/favicon.ico"]
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
