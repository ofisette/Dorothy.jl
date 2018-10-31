# Common XDR routines (for the XTC and TRR formats)

module XDR

export
		read_xdr_num, write_xdr_num, read_xdr_int, write_xdr_int,
		read_xdr_float, size_of_xdr_int, size_of_xdr_ints!, decode_xdr_bits!,
		decode_xdr_ints!

@inline read_xdr_num(io::IO, T::Type) = ntoh(read(io, T))

@inline write_xdr_num(io::IO, x::Number, T::Type) =
		write(io, hton(T(x)))

@inline read_xdr_int(io::IO) = read_xdr_num(io, Int32)

@inline write_xdr_int(io::IO, x::Integer) = write_xdr_num(io, x, Int32)

@inline read_xdr_float(io::IO) = read_xdr_num(io, Float32)

function size_of_xdr_int(size::T) where {T<:Integer}
	local num::T, num_of_bits::T
    num = 1
    num_of_bits = 0
    while (size >= num && num_of_bits < 32)
		num_of_bits += 1
		num <<= 1
	end
    num_of_bits
end

function size_of_xdr_ints!(bytes::Vector{T}, sizes::Vector{T}) where
		{T<:Integer}
	local num_of_bytes::T, num_of_bits::T, tmp::T, bytecnt::T, num::T,
			last_byte::T
    num_of_bytes = 1
    num_of_bits = 0
    bytes[1] = 1
	for size::Int in sizes
		tmp = 0
		bytecnt = 1
		while bytecnt <= num_of_bytes
			tmp = bytes[bytecnt] * size + tmp
			bytes[bytecnt] = tmp & 0xff
			tmp >>= 8
			bytecnt += 1
		end
		while tmp != 0
			bytes[bytecnt] = tmp & 0xff
			bytecnt += 1
			tmp >>= 8
		end
		num_of_bytes = bytecnt - 1
	end
    num = 1
	last_byte = bytes[num_of_bytes]
	num_of_bytes -= 1
    while last_byte >= num
		num_of_bits += 1
		num *= 2
	end
    num_of_bits + num_of_bytes * 8
end

function decode_xdr_bits!(state::Vector{T}, buf::Vector{UInt8},
		num_of_bits::T) where {T<:Integer}
	local mask::T, cnt::T, lastbits::T, lastbyte::T, num::T

	mask = (1 << num_of_bits) - 1
	cnt = state[1]
	lastbits = state[2]
	lastbyte = state[3]
	num = 0

	while num_of_bits >= 8
		lastbyte = (lastbyte << 8) | buf[cnt+1]
		cnt += 1
		num |= (lastbyte >> lastbits) << (num_of_bits - 8)
		num_of_bits -= 8
	end
	if num_of_bits > 0
		if lastbits < num_of_bits
			lastbits += 8
			lastbyte = (lastbyte << 8) | buf[cnt+1]
			cnt += 1
		end
		lastbits -= num_of_bits
		num |= (lastbyte >> lastbits) & ((1 << num_of_bits) - 1)
	end
	num &= mask
	state[1] = cnt
	state[2] = lastbits
	state[3] = lastbyte
	num
end

function decode_xdr_ints!(nums::Vector{T}, bytes::Vector{T},
		state::Vector{T}, buf::Vector{UInt8}, num_of_ints::T, num_of_bits::T,
		sizes::Vector{T}, pos::T) where {T<:Integer}
	@inbounds begin
		local num_of_bytes::T, i::T, j::T, p::T
		bytes[2] = 0
		bytes[3] = 0
		bytes[4] = 0
		num_of_bytes = 0
		while num_of_bits > 8
			bytes[num_of_bytes+1] = decode_xdr_bits!(state, buf, T(8))
			num_of_bytes += 1
			num_of_bits -= T(8)
		end
		if num_of_bits > 0
			bytes[num_of_bytes+1] = decode_xdr_bits!(state, buf, num_of_bits)
			num_of_bytes += 1
		end
		i = num_of_ints - 1
		while i > 0
			num::Int = 0
			j = num_of_bytes - 1
			while j >= 0
				num = (num << 8) | bytes[j+1]
				p = unsafe_trunc(Int, num / sizes[i+1])
				bytes[j+1] = p
				num -= p * sizes[i+1]
				j -= 1
			end
			nums[i+pos] = num
			i -= 1
		end
		nums[pos] = bytes[1] | (bytes[2] << 8) | (bytes[3] << 16) |
				(bytes[4] << 24)
	end
	nothing
end

end # module
