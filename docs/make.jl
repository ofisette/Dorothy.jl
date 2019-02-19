using Documenter
using Dorothy

docargs = (
	sitename="Dorothy.jl",
	pages = [
		"Home" => "index.md",
		"Introduction" => "intro.md",
		"Manual" => [],
		"Development" => [],
		"Appendices" => [
			"References" => "refs.md"
			"Software license" => "license.md"
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
