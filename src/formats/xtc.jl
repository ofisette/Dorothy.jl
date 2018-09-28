# XTC molecular trajectory file format

module XTC

using ..Dorothy
using ..Dorothy.XDR
using Formats
using FormatStreams

export register_xtc, readxtc!

const XTCInt = Int
const xtc_format_code = 1995
const xtc_magicints =
		XTCInt[0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
		80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
		1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
		16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031,
		131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
		832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
		4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216]
const xtc_firstidx = XTCInt(9)
const xtc_lastidx = XTCInt(length(xtc_magicints))

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
	uic::Vector{XTCInt}
	fic::Vector{XTCInt}
	bytes::Vector{XTCInt}
	minint::Vector{XTCInt}
	maxint::Vector{XTCInt}
	sizeint::Vector{XTCInt}
	bitsizeint::Vector{XTCInt}
	sizesmall::Vector{XTCInt}
	prevcoord::Vector{XTCInt}
	state::Vector{XTCInt}
	cell::Matrix{Float64}

	function XTCBuffer(nparticles::Integer)
		cic = Vector{UInt8}(undef, floor(Int, 3*nparticles*1.5))
		uic = Vector{UInt8}(undef, 3*nparticles)
		fic = Vector{UInt8}(undef, 3*nparticles)
		bytes = Vector{XTCInt}(undef, 32)
		minint = Vector{XTCInt}(undef, 3)
		maxint = Vector{XTCInt}(undef, 3)
		sizeint = Vector{XTCInt}(undef, 3)
		bitsizeint = Vector{XTCInt}(undef, 3)
		sizesmall = Vector{XTCInt}(undef, 3)
		prevcoord = Vector{XTCInt}(undef, 3)
		state = Vector{XTCInt}(undef, 3)
		cell = Matrix{Float32}(undef, 3,3)
		new(cic, uic, fic, bytes, minint, maxint, sizeint, bitsizeint,
				sizesmall, prevcoord, state, cell)
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
	model.header.cell = buffer.cell
	skip(io, 4)
	R = get!(model, :R, undef)
	if meta.nparticles > 9
		read3dfcoord!(io, R, meta, buffer)
	else
		readxdrcoord!(io, R, meta.nparticles)
	end
	model
end

function read3dfcoord!(io::IO, R::AbstractArray{<:Real}, meta::XTCMetadata,
		buffer::XTCBuffer)
	local nparticles::XTCInt, smallidx::XTCInt, nbytes::XTCInt, bitsize::XTCInt,
			tmp::XTCInt, maxidx::XTCInt, minidx::XTCInt, smaller::XTCInt,
			smallnum::XTCInt, larger::XTCInt, remainder::XTCInt, i::XTCInt,
			j::XTCInt, k::XTCInt, run::XTCInt, d::XTCInt, x::XTCInt, y::XTCInt,
			z::XTCInt, flag::XTCInt, is_smaller::XTCInt, magicints::Vector{XTCInt},
			firstidx::XTCInt, lastidx::XTCInt

	magicints = xtc_magicints
	firstidx = xtc_firstidx
	lastidx = xtc_lastidx

	nparticles = meta.nparticles
	precision = meta.precision
	cic = buffer.cic
	uic = buffer.uic
	fic = buffer.fic
	bytes = buffer.bytes
	minint = buffer.minint
	maxint = buffer.maxint
	sizeint = buffer.sizeint
	bitsizeint = buffer.bitsizeint
	sizesmall = buffer.sizesmall
	prevcoord = buffer.prevcoord
	state = buffer.state

	skip(io, 4)
	minint[1] = read_xdr_int(io)
	minint[2] = read_xdr_int(io)
	minint[3] = read_xdr_int(io)
	maxint[1] = read_xdr_int(io)
	maxint[2] = read_xdr_int(io)
	maxint[3] = read_xdr_int(io)
	smallidx = read_xdr_int(io)
	nbytes = read_xdr_int(io)
	if smallidx == 0
		error("not sure what this error is...")
	end
	if nbytes == 0
		error("zero-length compressed coordinates")
	end

	sizeint[1] = maxint[1] - minint[1] + 1
	sizeint[2] = maxint[2] - minint[2] + 1
	sizeint[3] = maxint[3] - minint[3] + 1
	if (sizeint[1] | sizeint[2] | sizeint[3]) > 0xffffff
		bitsizeint[1] = size_of_xdr_int(sizeint[1])
		bitsizeint[2] = size_of_xdr_int(sizeint[2])
		bitsizeint[3] = size_of_xdr_int(sizeint[3])
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
	sizesmall[1] = magicints[smallidx+1]
	sizesmall[2] = magicints[smallidx+1]
	sizesmall[3] = magicints[smallidx+1]
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
			decode_xdr_ints!(uic, bytes, state, cic, XTCInt(3), bitsize,
			sizeint, x)
		end
		i += 1

		uic[x] += minint[1]
		uic[y] += minint[2]
		uic[z] += minint[3]

		prevcoord[1] = uic[x]
		prevcoord[2] = uic[y]
		prevcoord[3] = uic[z]

		flag = decode_xdr_bits!(state, cic, XTCInt(1))
		is_smaller = 0
		if flag == 1
			run = decode_xdr_bits!(state, cic, XTCInt(5))
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
				decode_xdr_ints!(uic, bytes, state, cic, XTCInt(3), smallidx,
						sizesmall, x)
				i += 1
				uic[x] += prevcoord[1] - smallnum
				uic[y] += prevcoord[2] - smallnum
				uic[z] += prevcoord[3] - smallnum
				if k == 0
					tmp = uic[x]
					uic[x] = prevcoord[1]
					prevcoord[1] = tmp
					tmp = uic[y]
					uic[y] = prevcoord[2]
					prevcoord[2] = tmp
					tmp = uic[z]
					uic[z] = prevcoord[3]
					prevcoord[3] = tmp
					fic[j] = prevcoord[1]
					j += 1
					fic[j] = prevcoord[2]
					j += 1
					fic[j] = prevcoord[3]
					j += 1
				else
					prevcoord[1] = uic[x]
					prevcoord[2] = uic[y]
					prevcoord[3] = uic[z]
				end
				fic[j] = prevcoord[1]
				j += 1
				fic[j] = prevcoord[2]
				j += 1
				fic[j] = prevcoord[3]
				j += 1
				k += 3
			end
		else
			fic[j] = prevcoord[1]
			j += 1
			fic[j] = prevcoord[2]
			j += 1
			fic[j] = prevcoord[3]
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
		sizesmall[1] = sizesmall[2] = sizesmall[3] = magicints[smallidx + 1]
	end

	for i in eachindex(R)
		R[i] = fic[i] * 10 / precision
	end
	R
end

function readxdrcoord!(io::IO, R::AbstractArray{<:Real}, nparticles::Integer)
	for i = 1:3*nparticles
		R[i] = read_xdr_float(io)
	end
	R
end

end # module
