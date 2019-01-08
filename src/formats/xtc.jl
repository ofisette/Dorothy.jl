# XTC molecular trajectory file format

module XTC

using ..Dorothy
using ..Dorothy.Geometry
using ..Dorothy.PBC
using ..Dorothy.XDR
using Formats
using FormatStreams
using StaticArrays

export register_xtc, readxtc!

const xtc_format_code = 1995
const xtc_magicints =
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
		80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
		1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
		16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031,
		131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
		832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
		4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216]
const xtc_firstidx = 9
const xtc_lastidx = length(xtc_magicints)

function register_xtc()
	Formats.addformat("trajectory/x-xtc")
	Formats.addextension("trajectory/x-xtc", ".xtc")
	Formats.addsignature("trajectory/x-xtc", [0x00,0x00,0x07,0xcb])
	Formats.addreader("trajectory/x-xtc", DorothyIO())
	Formats.addwriter("trajectory/x-xtc", DorothyIO())
	FormatStreams.addstreamer("trajectory/x-xtc", DorothyIO())
end

FormatStreams.streamf(::DorothyIO, ::MIME"trajectory/x-xtc", io::T, args...;
		kwargs...) where {T<:IO} = XTCStream{T}(io, args...; kwargs...)

Base.read(::DorothyIO, ::MIME"trajectory/x-xtc", io::IO, args...; kwargs...) =
		readxtc!(io, MolecularModel(), args...; kwargs...)

Base.read!(::DorothyIO, ::MIME"trajectory/x-xtc", io::IO, model::MolecularModel,
		args...; kwargs...) = readxtc!(io, model, args...; kwargs...)

Base.write(::DorothyIO, ::MIME"trajectory/x-xtc", io::IO, model::MolecularModel,
		args...; kwargs...) = writextc(io, model, args...; kwargs...)

struct XTCMetadata
	nparticles::Int
	precision::Float64

	function XTCMetadata(nparticles::Integer, precision::Real = 1000.0)
		@boundscheck begin
			nparticles >= 0 || error("expected positive number of particles")
			precision >= 0.0 || error("expected positive precision")
		end
		if nparticles <= 9
			precision = 0.0
		end
		new(nparticles, precision)
	end
end

function read_xtc_metadata(io::IO)
	skip(io, 4)
	nparticles = read_xdr_int(io)
	if nparticles <= 9
		precision = 0.0
	else
		skip(io, 48)
		precision = read_xdr_float(io)
	end
	XTCMetadata(nparticles, precision)
end

function peek_xtc_metadata(io::IO)
	mark(io)
	meta = read_xtc_metadata(io)
	reset(io)
	meta
end

struct XTCBuffer
	cic::Vector{UInt8}
	uic::Vector{Int}
	fic::Vector{Int}
	bytes::MVector{32,Int}
	state::MVector{3,Int}
	cell::MMatrix{3,3,Float64}

	function XTCBuffer(nparticles::Integer)
		cic = Vector{UInt8}(undef, floor(Int, 3*nparticles*1.5))
		uic = Vector{UInt8}(undef, 3*nparticles)
		fic = Vector{UInt8}(undef, 3*nparticles)
		bytes = MVector{32,Int}(undef)
		state = MVector{3,Int}(undef)
		cell = MMatrix{3,3,Float64}(undef)
		new(cic, uic, fic, bytes, state, cell)
	end
end

mutable struct XTCStream{T<:IO} <: MolecularTrajectory
	io::T
	frameindices::Vector{Int}
	pos::Int
	meta::XTCMetadata
	buffer::XTCBuffer

	function XTCStream{T}(io::T; nparticles::Integer = -1,
			precision::Real = 1000.0) where {T<:IO}
		origin = position(io)
		if eof(io)
			meta = XTCMetadata(nparticles, precision)
		else
			meta = peek_xtc_metadata(io)
		end
		new(io, [origin], 0, meta, XTCBuffer(meta.nparticles))
	end
end

Base.close(s::XTCStream) = close(s.io)

Base.eof(s::XTCStream) = eof(s.io)

Base.position(s::XTCStream) = s.pos

function Base.seek(s::XTCStream, pos::Integer)
	i = pos + 1
	@boundscheck begin
		i >= 0 || error("expected positive position")
		if i > length(s.frameindices)
			error("cannot seek to unknown position")
		end
	end
	seek(s.io, s.frameindices[i])
	s.pos = pos
	s
end

function Base.seekstart(s::XTCStream)
	seek(s.io, s.frameindices[1])
	s.pos = 0
	s
end

Base.read(s::XTCStream) = read!(s, MolecularModel())

function Base.read!(s::XTCStream, model::MolecularModel)
	readxtc!(s.io, model, s.meta, s.buffer)
	s.pos += 1
	if s.pos + 1 > length(s.frameindices)
		push!(s.frameindices, position(s.io))
	end
	model
end

#=
function Base.write(s::XTCStream, model::MolecularModel)
	@boundscheck eof(s.io) || error("writing must happen at end of trajectory")
	writextc(s.io, model, s.meta, s.buffer)
	s.pos += 1
	push!(s.frameindices, position(s.io))
	s.nframes += 1
end
=#

function Base.truncate(s::XTCStream, n::Integer)
	@boundscheck begin
		n > 0 || error("expected positive value")
		n + 1 > length(s.frameindices) &&
				error("cannot truncate to unknown position")
	end
	truncate(s.io, s.frameindices[n+1])
	deleteat!(s.frameindices, n+2:lastindex(s.frameindices))
	if s.pos > n
		s.pos = n
	end
	s
end

function readxtc!(io::IO, model::MolecularModel,
		meta::XTCMetadata = peek_xtc_metadata(io),
		buffer::XTCBuffer = XTCBuffer(meta.nparticles))
	skip(io, 8)
	resize!(model, meta.nparticles)
	model.header.step = read_xdr_int(io)
	model.header.time = read_xdr_float(io)
	for i = 1:9
		buffer.cell[i] = read_xdr_float(io) * 10.0
	end
	model.header.cell = pbccell(buffer.cell)
	skip(io, 4)
	R = get!(model, :R, undef)
	if meta.nparticles > 9
		read3dfcoord!(io, R, meta, buffer)
	else
		readxdrcoord!(io, R, meta.nparticles)
	end
	model
end

function read3dfcoord!(io::IO, R::AbstractVector{Vector3D}, meta::XTCMetadata,
		buffer::XTCBuffer)

	magicints = xtc_magicints
	firstidx = xtc_firstidx
	lastidx = xtc_lastidx

	nparticles = meta.nparticles
	precision = meta.precision
	cic = buffer.cic
	uic = buffer.uic
	fic = buffer.fic
	bytes = buffer.bytes
	state = buffer.state

	skip(io, 4)
	minint = SVector{3,Int}(read_xdr_int(io), read_xdr_int(io),
			read_xdr_int(io))
	maxint = SVector{3,Int}(read_xdr_int(io), read_xdr_int(io),
			read_xdr_int(io))
	smallidx = Int(read_xdr_int(io))
	nbytes = Int(read_xdr_int(io))
	if smallidx == 0
		error("not sure what this error is...")
	end
	if nbytes == 0
		error("zero-length compressed coordinates")
	end

	sizeint = maxint - minint .+ 1
	if (sizeint[1] | sizeint[2] | sizeint[3]) > 0xffffff
		bitsizeint = size_of_xdr_int.(sizeint)
		bitsize = 0
	else
		bitsize = size_of_xdr_ints!(bytes, sizeint)
	end

	tmp = smallidx + 8
	maxidx = lastidx < tmp ? lastidx : tmp
	minidx = maxidx - 8
	tmp = smallidx - 1
	tmp = firstidx > tmp ? firstidx : tmp
	smaller = div(magicints[tmp+1], 2)
	smallnum = div(magicints[smallidx+1], 2)
	sizesmall = SVector(magicints[smallidx+1], magicints[smallidx+1],
			magicints[smallidx+1])
	larger = magicints[maxidx+1]

	resize!(cic, nbytes)
	read!(io, cic)
	remainder = nbytes % 4
	if remainder != 0
		skip(io, 4 - remainder)
	end
	for i in 1:3
		state[i] = 0
	end

	i = 0
	j = 1
	k = 0
	run = 0
	while i < nparticles
		d = 3 * i
		x = d + 1
		y = d + 2
		z = d + 3
		if bitsize == 0
			uic[x] = decode_xdr_bits!(state, cic, bitsizeint[1])
			uic[y] = decode_xdr_bits!(state, cic, bitsizeint[2])
			uic[z] = decode_xdr_bits!(state, cic, bitsizeint[3])
		else
			decode_xdr_ints!(uic, bytes, state, cic, 3, bitsize, sizeint, x)
		end
		i += 1

		uic[x] += minint[1]
		uic[y] += minint[2]
		uic[z] += minint[3]

		prevcoordx = uic[x]
		prevcoordy = uic[y]
		prevcoordz = uic[z]

		flag = decode_xdr_bits!(state, cic, 1)
		is_smaller = 0
		if flag == 1
			run = decode_xdr_bits!(state, cic, 5)
			is_smaller = run % 3
			run -= is_smaller
			is_smaller -= 1
		end
		if run > 0
			d += 3
			x = d + 1
			y = d + 2
			z = d + 3
			k = 0
			while k < run
				decode_xdr_ints!(uic, bytes, state, cic, 3, smallidx,
						sizesmall, x)
				i += 1
				uic[x] += prevcoordx - smallnum
				uic[y] += prevcoordy - smallnum
				uic[z] += prevcoordz - smallnum
				if k == 0
					tmp = uic[x]
					uic[x] = prevcoordx
					prevcoordx = tmp
					tmp = uic[y]
					uic[y] = prevcoordy
					prevcoordy = tmp
					tmp = uic[z]
					uic[z] = prevcoordz
					prevcoordz = tmp
					fic[j] = prevcoordx
					j += 1
					fic[j] = prevcoordy
					j += 1
					fic[j] = prevcoordz
					j += 1
				else
					prevcoordx = uic[x]
					prevcoordy = uic[y]
					prevcoordz = uic[z]
				end
				fic[j] = uic[x]
				j += 1
				fic[j] = uic[y]
				j += 1
				fic[j] = uic[z]
				j += 1
				k += 3
			end
		else
			fic[j] = uic[x]
			j += 1
			fic[j] = uic[y]
			j += 1
			fic[j] = uic[z]
			j += 1
		end
		smallidx += is_smaller
		if is_smaller < 0
			smallnum = smaller
			if smallidx > firstidx
				smaller = div(magicints[smallidx], 2)
			else
				smaller = 0
			end
		elseif is_smaller > 0
			smaller = smallnum
			smallnum = div(magicints[smallidx+1], 2)
		end
		sizesmall = SVector(magicints[smallidx+1], magicints[smallidx+1],
				magicints[smallidx+1])
	end

	for i in eachindex(R)
		d = (i-1) * 3
		x = (fic[d+1] * 10) / precision
		y = (fic[d+2] * 10) / precision
		z = (fic[d+3] * 10) / precision
		R[i] = Vector3D(x,y,z)
	end
	R
end

function readxdrcoord!(io::IO, R::AbstractArray{<:Real}, nparticles::Integer)
	for i = 1:nparticles
		x = read_xdr_float(io)
		y = read_xdr_float(io)
		z = read_xdr_float(io)
		R[i] = Vector3D(x,y,z)
	end
	R
end

end # module
