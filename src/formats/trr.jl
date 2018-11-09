# TRR molecular trajectory file format

module TRR

using ..Dorothy
using ..Dorothy.Geometry
using ..Dorothy.XDR
using Formats
using FormatStreams

export register_trr, readtrr!, writetrr

const trr_format_signature = "\0\0\a\xc9\0\0\0\r\0\0\0\fGMX_trn_file"
const trr_precision_types = [Float32, Float64]

function register_trr()
	Formats.addformat("trajectory/x-trr")
	Formats.addextension("trajectory/x-trr", ".trr")
	Formats.addsignature("trajectory/x-trr", trr_format_signature)
	Formats.addreader("trajectory/x-trr", DorothyIO())
	Formats.addwriter("trajectory/x-trr", DorothyIO())
	FormatStreams.addstreamer("trajectory/x-trr", DorothyIO())
end

FormatStreams.streamf(::DorothyIO, ::MIME"trajectory/x-trr", io::T, args...;
		kwargs...) where {T<:IO} = TRRStream{T}(io, args...; kwargs...)

Base.read(::DorothyIO, ::MIME"trajectory/x-trr", io::IO, args...; kwargs...) =
		readtrr!(io, MolecularModel(), args...; kwargs...)

Base.read!(::DorothyIO, ::MIME"trajectory/x-trr", io::IO, model::MolecularModel,
		args...; kwargs...) = readtrr!(io, model, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"trajectory/x-trr", io::IO,
		model::ParticleCollection, args...; kwargs...) =
		writetrr(io, model, args...; kwargs...)

struct TRRMetadata
	irsize::Int
	esize::Int
	cellsize::Int
	virialsize::Int
	pressuresize::Int
	topsize::Int
	symsize::Int
	possize::Int
	Vsize::Int
	Fsize::Int
	nparticles::Int
	precision::Type
	framesize::Int

	function TRRMetadata(irsize::Integer, esize::Integer, cellsize::Integer,
			virialsize::Integer, pressuresize::Integer, topsize::Integer,
			symsize::Integer, possize::Integer, Vsize::Integer,
			Fsize::Integer, nparticles::Integer, precision::Type)
		framesize = 76 + 2 * sizeof(precision) + cellsize + virialsize +
				pressuresize + possize + Vsize + Fsize
		new(irsize, esize, cellsize, virialsize, pressuresize, topsize,
				symsize, possize, Vsize, Fsize, nparticles, precision,
				framesize)
	end
end

function read_trr_metadata(io::IO)
	skip(io, 24)
	irsize = read_xdr_int(io)
	esize = read_xdr_int(io)
	cellsize = read_xdr_int(io)
	virialsize = read_xdr_int(io)
	pressuresize = read_xdr_int(io)
	topsize = read_xdr_int(io)
	symsize = read_xdr_int(io)
	possize = read_xdr_int(io)
	Vsize = read_xdr_int(io)
	Fsize = read_xdr_int(io)
	nparticles = read_xdr_int(io)
	precision = guess_trr_precision([cellsize, virialsize, pressuresize],
			[possize, Vsize, Fsize], nparticles)
	TRRMetadata(irsize, esize, cellsize, virialsize, pressuresize,
			topsize, symsize, possize, Vsize, Fsize, nparticles,
			precision)
end

function peek_trr_metadata(io::IO)
	mark(io)
	meta = read_trr_metadata(io)
	reset(io)
	meta
end

skip_trr_metadata(io::IO) = skip(io, 68)

function write_trr_metadata(io::IO, meta::TRRMetadata)
	write(io, trr_format_signature)
	write_xdr_int(io, meta.irsize)
	write_xdr_int(io, meta.esize)
	write_xdr_int(io, meta.cellsize)
	write_xdr_int(io, meta.virialsize)
	write_xdr_int(io, meta.topsize)
	write_xdr_int(io, meta.symsize)
	write_xdr_int(io, meta.pressuresize)
	write_xdr_int(io, meta.possize)
	write_xdr_int(io, meta.Vsize)
	write_xdr_int(io, meta.Fsize)
	write_xdr_int(io, meta.nparticles)
end

function trr_metadata_for_model(model::MolecularModel)
	properties = Symbol[]
	for key in [:cell, :virial, :pressure, :R, :V, :F]
		if haskey(model, key)
			push!(properties, key)
		end
	end
	trr_metadata_from_vals(length(model), Float32, properties)
end

function trr_metadata_from_vals(nparticles::Integer, precision::Type = Float32,
		properties::AbstractVector{Symbol} = [:cell, :R])
	@boundscheck begin
		nparticles >= 0 || error("expected positive number of particles")
		precision in trr_precision_types ||
				error("invalid TRR trajectory precision")
	end
	cellsize = (:cell in properties ? 3*3 * sizeof(precision) : 0)
	virialsize = (:virial in properties ? 3*3 * sizeof(precision) : 0)
	pressuresize = (:pressure in properties ? 3*3 * sizeof(precision) : 0)
	possize = (:R in properties ? 3*nparticles * sizeof(precision) : 0)
	Vsize = (:V in properties ? 3*nparticles * sizeof(precision) : 0)
	Fsize = (:F in properties ? 3*nparticles * sizeof(precision) : 0)
	TRRMetadata(0, 0, cellsize, virialsize, pressuresize, 0, 0, possize,
			Vsize, Fsize, nparticles, precision)
end

function guess_trr_precision(header_array_sizes::AbstractVector{<:Integer},
		body_array_sizes::AbstractVector{<:Integer}, nparticles::Integer)
	for array_size in header_array_sizes
		for T in trr_precision_types
			if array_size == 3 * 3 * sizeof(T)
				return T
			end
		end
	end
	for array_size in body_array_sizes
		for T in trr_precision_types
			if array_size == 3 * nparticles * sizeof(T)
				return T
			end
		end
	end
	error("could not determine TRR trajectory precision")
end

struct TRRBuffer{T<:Real}
	cell::Matrix{T}
	cellf::Matrix{Float64}
	virial::Matrix{T}
	virialf::Matrix{Float64}
	pressure::Matrix{T}
	pressuref::Matrix{Float64}
	R::Matrix{T}
	V::Matrix{T}
	F::Matrix{T}

	function TRRBuffer{T}(meta::TRRMetadata) where {T<:Real}
		n = meta.cellsize > 0 ? 3 : 0
		cell = Matrix{T}(undef, 3,n)
		cellf = Matrix{Float64}(undef, 3,n)
		n = meta.virialsize > 0 ? 3 : 0
		virial = Matrix{T}(undef, 3,n)
		virialf = Matrix{Float64}(undef, 3,n)
		n = meta.pressuresize > 0 ? 3 : 0
		pressure = Matrix{T}(undef, 3,n)
		pressuref = Matrix{Float64}(undef, 3,n)
		n = meta.possize > 0 ? meta.nparticles : 0
		R = Matrix{T}(undef, 3,n)
		n = meta.Vsize > 0 ? meta.nparticles : 0
		V = Matrix{T}(undef, 3,n)
		n = meta.Fsize > 0 ? meta.nparticles : 0
		F = Matrix{T}(undef, 3,n)
		new(cell, cellf, virial, virialf, pressure, pressuref, R, V, F)
	end
end

TRRBuffer(meta::TRRMetadata) = TRRBuffer{meta.precision}(meta)

mutable struct TRRStream{T<:IO} <: MolecularTrajectory
	io::T
	origin::Int
	pos::Int
	nframes::Int
	meta::TRRMetadata
	buffer::TRRBuffer

	function TRRStream{T}(io::T; nparticles::Integer = -1,
			precision::Type = Float32,
			properties::AbstractVector{Symbol} = [:cell, :R]) where {T<:IO}
		origin = position(io)
		if eof(io)
			nparticles >= 0 || error("expected positive number of particles")
			meta = trr_metadata_from_vals(nparticles, precision, properties)
		else
			meta = read_trr_metadata(io)
		end
		seekend(io)
		nframes = (position(io) - origin) ÷ meta.framesize
		seek(io, origin)
		new(io, origin, 0, nframes, meta, TRRBuffer(meta))
	end
end

Base.close(s::TRRStream) = close(s.io)

Base.length(s::TRRStream) = s.nframes

Base.eof(s::TRRStream) = (s.pos >= s.nframes)

Base.position(s::TRRStream) = s.pos

function Base.seek(s::TRRStream, pos::Integer)
	@boundscheck checkbounds(0:s.nframes, pos)
	seek(s.io, s.origin + s.meta.framesize * pos)
	s.pos = pos
	s
end

function Base.seekstart(s::TRRStream)
	seek(s.io, s.origin)
	s.pos = 0
	s
end

function Base.seekend(s::TRRStream)
	seekend(s.io)
	s.pos = s.nframes
	s
end

Base.read(s::TRRStream) = read!(s, MolecularModel())

function Base.read!(s::TRRStream, model::MolecularModel)
	readtrr!(s.io, model, s.meta, s.buffer)
	s.pos += 1
	model
end

function Base.write(s::TRRStream, model::ParticleCollection)
	pos0 = position(s.io)
	writetrr(s.io, model, s.meta, s.buffer)
	s.pos += 1
	s.nframes += 1
	position(s.io) - pos0
end

function Base.truncate(s::TRRStream, n::Integer)
	@boundscheck n >= 0 || error("expected positive value")
	truncate(s.io, n * s.meta.framesize)
	s.nframes = n
	if s.pos > n
		s.pos = n
	end
	s
end

function readtrr!(io::IO, model::MolecularModel,
		meta::TRRMetadata = peek_trr_metadata(io),
		buffer::TRRBuffer = TRRBuffer(meta))
	skip_trr_metadata(io)
	resize!(model, meta.nparticles)
	model.header.step = read_xdr_int(io)
	skip(io, 4)
	model.header.time = read_xdr_num(io, meta.precision)
	model.header.lambda = read_xdr_num(io, meta.precision)
	if meta.cellsize != 0
		read!(io, buffer.cell)
		decode_trr_array!(buffer.cellf, buffer.cell, 10.0)
		model.header.cell = buffer.cellf
	end
	if meta.virialsize != 0
		read!(io, buffer.virial)
		decode_trr_array!(buffer.virialf, buffer.virial)
		model.header.virial = buffer.virialf
	end
	if meta.pressuresize != 0
		read!(io, buffer.pressure)
		decode_trr_array!(buffer.pressuref, buffer.pressure, 1.0/1000.0)
		model.header.pressure = buffer.pressuref
	end
	if meta.possize != 0
		read!(io, buffer.R)
		decode_trr_array!(get!(model, :R, undef), buffer.R, 10.0)
	end
	if meta.Vsize != 0
		read!(io, buffer.V)
		decode_trr_array!(get!(model, :V, undef), buffer.V, 10.0)
	end
	if meta.Fsize != 0
		read!(io, buffer.F)
		decode_trr_array!(get!(model, :F, undef), buffer.F, 1/10.0)
	end
	model
end

function decode_trr_array!(dest::AbstractArray{<:Real},
		src::AbstractArray{<:Real}, factor::Real = 1.0)
	@inbounds @fastmath @simd for i in eachindex(dest)
		dest[i] = ntoh(src[i]) * factor
	end
	dest
end

function writetrr(io::IO, model::ParticleCollection,
		meta::TRRMetadata = trr_metadata_for_model(model),
		buffer::TRRBuffer = TRRBuffer(meta))
	write_trr_metadata(io, meta)
	write_xdr_int(io, model.header.step)
	write_xdr_int(io, 0)
	write_xdr_num(io, model.header.time, meta.precision)
	write_xdr_num(io, model.header.lambda, meta.precision)
	if meta.cellsize != 0
		encode_trr_array!(buffer.cell, model.header.cell, 10.0)
		write(io, buffer.cell)
	end
	if meta.virialsize != 0
		encode_trr_array!(buffer.virial, model.header.virial)
		write(io, buffer.virial)
	end
	if meta.pressuresize != 0
		encode_trr_array!(buffer.pressure, model.header.pressure, 1.0/1000.0)
		write(io, buffer.pressure)
	end
	if meta.possize != 0
		encode_trr_array!(buffer.R, model.R, 10.0)
		write(io, buffer.R)
	end
	if meta.Vsize != 0
		encode_trr_array!(buffer.V, model.V, 10.0)
		write(io, buffer.V)
	end
	if meta.Fsize != 0
		encode_trr_array!(buffer.F, model.F, 1/10.0)
		write(io, buffer.F)
	end
end

function encode_trr_array!(dest::AbstractArray{T},
		src::AbstractArray{<:Real}, factor::Real = 1.0) where {T<:Real}
	@inbounds @fastmath @simd for i in eachindex(dest)
		dest[i] = hton(T(src[i] / factor))
	end
	dest
end

end # module
