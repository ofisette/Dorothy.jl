# Algorithms on basic atomic properties

module Atoms

using ..Dorothy
using ..Dorothy.Properties
using ..Dorothy.Utils

export
		inferelement, inferelements!, infermissingelements,
		infermissingelements!, infermass, infermasses!

function inferelement(name::AbstractString, resname::AbstractString)
	if ismonatomicion(resname)
		titlecase(resname)
	elseif isempty(name)
		error("could not infer element from name or resname")
	else
		name[1:1]
	end
end

inferelements!(model::ParticleCollection) =
		inferelements!(get!(model, :elements, undef), model)

function inferelements!(elements::AbstractVector{<:AbstractString},
		model::ParticleCollection)
	@boundscheck length(elements) == length(model) ||
			error("size mismatch between model and output array")
	inferelements!(elements, model.names, model.resnames)
end

function inferelements!(elements::AbstractVector{<:AbstractString},
		names::AbstractVector{<:AbstractString},
		resnames::AbstractVector{<:AbstractString})
	@boundscheck begin
		length(elements) == length(names) == length(resnames) ||
				error("size mismatch between property arrays")
	end
	for i in eachindex(elements)
		elements[i] = inferelement(names[i], resnames[i])
	end
	elements
end

function infermissingelements(model::ParticleCollection)
	elements = get(model, :elements) do
		Repeated("", length(model))
	end
	infermissingelements!(collect(elements), model)
end

infermissingelements!(model::ParticleCollection) =
		infermissingelements!(model.elements, model)

function infermissingelements!(elements::AbstractVector{<:AbstractString},
		model::ParticleCollection)
	@boundscheck length(elements) == length(model) ||
			error("size mismatch between model and output array")
	inferelements!(elements, model.names, model.resnames)
end

function infermissingelements!(elements::AbstractVector{<:AbstractString},
		names::AbstractVector{<:AbstractString},
		resnames::AbstractVector{<:AbstractString})
	@boundscheck begin
		length(elements) == length(names) == length(resnames) ||
				error("size mismatch between property arrays")
	end
	for i in eachindex(elements)
		if elements[i] == ""
			elements[i] = inferelement(names[i], resnames[i])
		end
	end
	elements
end

function infermass(name::AbstractString, element::AbstractString)
	get(standard_atomic_weights, element) do
		if isvsite(name)
			0.0
		else
			error("could not infer mass from name or element")
		end
	end
end

infermasses!(model::ParticleCollection) =
		infermasses!(get!(model, :masses, undef), model)

function infermasses!(masses::AbstractVector{<:Real}, model::ParticleCollection)
	@boundscheck length(masses) == length(model) ||
			error("size mismatch between model and output array")
	infermasses!(masses, model.names,
			get(model, :elements, inferelements!(model)))
end

function infermasses!(masses::AbstractVector{<:Real},
		names::AbstractVector{<:AbstractString},
		elements::AbstractVector{<:AbstractString})
	@boundscheck begin
		length(masses) == length(names) == length(elements) ||
				error("size mismatch between property arrays")
	end
	for i in eachindex(masses)
		masses[i] = infermass(names[i], elements[i])
	end
	masses
end

end # module
