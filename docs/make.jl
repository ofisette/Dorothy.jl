using Documenter
using Dorothy

makedocs(
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
	]
)
