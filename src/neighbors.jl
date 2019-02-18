module Neighbors

using ..Dorothy.Geometry
using ..Dorothy.PBC

export
		NeighborGridBuffer, NeighborGrid, findnear, findnear!, findwithin,
		findwithin!, NeighborDistance, neighborlist, neighborlist!,
		NeighborTableBuffer, NeighborTable

mutable struct NeighborGridBuffer
	cells::Array{Vector{Int},3}
	I::Vector{Tuple{Int,Int,Int}}
end

NeighborGridBuffer() = NeighborGridBuffer(Array{Vector{Int},3}(undef, (0,0,0)),
		Tuple{Int,Int,Int}[])

function prepare!(buffer::NeighborGridBuffer, N::Tuple{Integer,Integer,Integer})
	nx, ny, nz = N
	mx, my, mz = size(buffer.cells)
	if nx > mx || ny > my || nz > mz
		buffer.cells = newcells = Array{Vector{Int},3}(undef, N .+ 1)
		for i in eachindex(newcells)
			newcells[i] = Int[]
		end
	else
		for oldcell in buffer.cells
			empty!(oldcell)
		end
	end
	empty!(buffer.I)
	buffer
end

abstract type NeighborGrid end

struct NonperiodicNeighborGrid <: NeighborGrid
	cells::Array{Vector{Int},3}
	I::Vector{Tuple{Int,Int,Int}}
	N::Tuple{Int,Int,Int}
	O::Vector3D
	D::Vector3D
end

struct PeriodicNeighborGrid <: NeighborGrid
	cells::Array{Vector{Int},3}
	I::Vector{Tuple{Int,Int,Int}}
	N::Tuple{Int,Int,Int}
	D::Vector3D
end

NeighborGrid(N::Tuple{Integer,Integer,Integer}, O::AbstractVector{<:Real},
		D::AbstractVector{<:Real}, buffer::NeighborGridBuffer) =
		NonperiodicNeighborGrid(buffer.cells, buffer.I, N, O, D)

NeighborGrid(N::Tuple{Integer,Integer,Integer}, D::AbstractVector{<:Real},
		buffer::NeighborGridBuffer) =
		PeriodicNeighborGrid(buffer.cells, buffer.I, N, D)

function NeighborGrid(Rw::AbstractVector{Vector3D},
		Kw::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, buffer = NeighborGridBuffer())
	N, = params = neighborgridparams(Rw, Kw, cell, dmax)
	prepare!(buffer, N)
	nbrgrid = NeighborGrid(params..., buffer)
	append!(nbrgrid, Kw)
end

function neighborgridparams(Rw::AbstractVector{Vector3D},
		Kw::AbstractVector{Vector3D}, ::Nothing, dmax::Real)
	D = Vector3D(dmax, dmax, dmax)
	bb = extent(Rw)
	O = minimum(bb) - 2*D
	Dbb = dims(bb)
	nx = ceil(Int, Dbb.x / D.x) + 4
	ny = ceil(Int, Dbb.y / D.y) + 4
	nz = ceil(Int, Dbb.z / D.z) + 4
	(nx,ny,nz), O, D
end

function neighborgridparams(Rw::AbstractVector{Vector3D},
		Kw::AbstractVector{Vector3D}, cell::TriclinicPBC, dmax::Real)
	D = neighborgriddims(cell, dmax)
	nx = max(1, floor(Int, 1 / D.x))
	ny = max(1, floor(Int, 1 / D.y))
	nz = max(1, floor(Int, 1 / D.z))
	D = Vector3D(1.0/nx, 1.0/ny, 1.0/nz)
	(nx,ny,nz), D
end

function neighborgriddims(cell::TriclinicPBC, dmax::Real)
	rmax = norm(inv(cell) * Vector3D(dmax, dmax, dmax))
	Vector3D(rmax, rmax, rmax)
end

function neighborgriddims(cell::RhombododecahedralPBC, dmax::Real)
	rmax = √(2*dmax^2)
	gridcell = RhombododecahedralCell(rmax)
	inv(cell) * dims(gridcell)
end

neighborgriddims(cell::OrthorhombicPBC, dmax::Real) =
		inv(cell) * Vector3D(dmax, dmax, dmax)

Base.length(nbrgrid::NeighborGrid) = length(nbrgrid.I)

Base.size(nbrgrid::NeighborGrid) = nbrgrid.N

function Base.sizehint!(nbrgrid::NeighborGrid, n::Integer,
		ncell::Integer = cld(n, prod(size(nbrgrid))))
	@boundscheck n >= 0 && ncell >= 0 || error("expected positive length")
	sizehint!(nbrgrid.I, n)
	for cell in nbrgrid.cells
		sizehint!(cell, ncell)
	end
	nbrgrid
end

findcell(nbrgrid::NeighborGrid, i::Integer) = nbrgrid.I[i]

function findcell(nbrgrid::NonperiodicNeighborGrid, R::Vector3D)
	x = floor(Int, (R.x - nbrgrid.O.x) / nbrgrid.D.x) + 1
	y = floor(Int, (R.y - nbrgrid.O.y) / nbrgrid.D.y) + 1
	z = floor(Int, (R.z - nbrgrid.O.z) / nbrgrid.D.z) + 1
	x, y, z
end

function findcell(nbrgrid::PeriodicNeighborGrid, Rw::Vector3D)
	x = floor(Int, Rw.x / nbrgrid.D.x) + 1
	y = floor(Int, Rw.y / nbrgrid.D.y) + 1
	z = floor(Int, Rw.z / nbrgrid.D.z) + 1
	x, y, z
end

function Base.push!(nbrgrid::NeighborGrid,
		(x,y,z)::Tuple{Integer,Integer,Integer})
	nx, ny, nz = size(nbrgrid)
	!(x in 1:nx) || !(y in 1:ny) || !(z in 1:nz) &&
			error("position is outside the grid")
	push!(nbrgrid.cells[x,y,z], length(nbrgrid.I) + 1)
	push!(nbrgrid.I, (x,y,z))
	nbrgrid
end

Base.push!(nbrgrid::NeighborGrid, R::Vector3D) =
		push!(nbrgrid, findcell(nbrgrid, R))

function Base.append!(nbrgrid::NeighborGrid, R::AbstractVector{Vector3D})
	sizehint!(nbrgrid, length(nbrgrid) + length(R))
	for Ri in R
		push!(nbrgrid, Ri)
	end
	nbrgrid
end

findnear(R::AbstractVector{<:Vector3D}, i::Integer, dmax::Real) =
		findnear!(Int[], R, i, dmax)

findnear(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findnear!(Int[], R, i, cell, dmax)

findnear(Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findnear!(Int[], Rw, Kw, i, cell, dmax)

findnear!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		i::Integer, dmax::Real) = findnear!(dest, R, i, nothing, dmax)

findnear!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findnear!(dest, pbcpos(R, cell)..., i, cell, dmax)

function findnear!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real)
	nbrgrid = NeighborGrid(Rw, Kw, cell, dmax)
	findnear!(dest, nbrgrid, i)
end

findnear(R::AbstractVector{<:Vector3D}, Rref::AbstractVector{<:Real},
		dmax::Real) = findnear!(Int[], R, Rref, dmax)

findnear(R::AbstractVector{<:Vector3D}, Rref::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findnear!(Int[], R, Rref, cell, dmax)

findnear(Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		Rrefw::AbstractVector{<:Real}, Krefw::AbstractVector{<:Vector3D},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findnear!(Int[], Rw, Kw, Rrefw, Krefw, cell, dmax)

findnear!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		Rref::AbstractVector{<:Real}, dmax::Real) =
		findnear!(dest, R, Rref, nothing, dmax)

findnear!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		Rref::AbstractVector{<:Real}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) =
		findnear!(dest, pbcpos(R, cell)..., pbcpos(Rref, cell)..., cell, dmax)

findnear!(dest::AbstractVector{<:Integer}, Rw::AbstractVector{<:Vector3D},
		Kw::AbstractVector{<:Vector3D}, Rrefw::AbstractVector{<:Real},
		Krefw::AbstractVector{<:Vector3D}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) =
		findnear!(dest, Rw, Kw, Vector3D(Rrefw), Vector3D(Krefw), cell, dmax)

function findnear!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real)
	nbrgrid = NeighborGrid(Rw, Kw, cell, dmax)
	findnear!(dest, nbrgrid, Krefw)
end

findnear!(dest::AbstractVector{<:Integer}, nbrgrid::NeighborGrid, i::Integer) =
		findnear!(dest, nbrgrid, findcell(nbrgrid, i))

findnear!(dest::AbstractVector{<:Integer}, nbrgrid::NeighborGrid, R::Vector3D) =
		findnear!(dest, nbrgrid, findcell(nbrgrid, R))

function findnear!(dest::AbstractVector{<:Integer},
		nbrgrid::NonperiodicNeighborGrid,
		(x,y,z)::Tuple{Integer,Integer,Integer})
	empty!(dest)
	nx, ny, nz = size(nbrgrid)
	for xi = x-1:x+1
		for yi = y-1:y+1
			for zi = z-1:z+1
				if 1 <= xi <= nx && 1 <= yi <= ny && 1 <= zi <= nz
					for i in nbrgrid.cells[xi,yi,zi]
						push!(dest, i)
					end
				end
			end
		end
	end
	dest
end

function findnear!(dest::AbstractVector{<:Integer},
		nbrgrid::PeriodicNeighborGrid, (x,y,z)::Tuple{Integer,Integer,Integer})
	empty!(dest)
	nx, ny, nz = size(nbrgrid)
	for xi = x-1:x+1
		for yi = y-1:y+1
			for zi = z-1:z+1
				if xi < 1
					xi = nx
				elseif xi > nx
					xi = 1
				end
				if yi < 1
					yi = ny
				elseif yi > ny
					yi = 1
				end
				if zi < 1
					zi = nz
				elseif zi > nz
					zi = 1
				end
				for i in nbrgrid.cells[xi,yi,zi]
					push!(dest, i)
				end
			end
		end
	end
	dest
end

findwithin(R::AbstractVector{Vector3D}, i::Integer, dmax::Real) =
		findwithin!(Int[], R, nothing, dmax)

findwithin(R::AbstractVector{Vector3D}, i::Integer,
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findwithin!(Int[], pbcpos(R, cell)..., i, cell, dmax)

findwithin(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findwithin!(Int[], Rw, Kw, i, cell, dmax)

findwithin!(dest::AbstractVector{<:Integer}, R::AbstractVector{Vector3D},
		i::Integer, dmax::Real) = findwithin!(dest, R, i, nothing, dmax)

findwithin!(dest::AbstractVector{<:Integer}, R::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findwithin!(dest, pbcpos(R, cell)..., i, cell, dmax)

function findwithin!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D}, i::Integer,
		cell::Union{TriclinicPBC,Nothing}, dmax::Real)
	dmax2 = dmax^2
	Rwi, Kwi = Rw[i], Kw[i]
	empty!(dest)
	for j in findnear(Rw, Kw, i, cell, dmax)
		if sqmindist(Rw[j], Kw[j], Rwi, Kwi, cell) <= dmax2
			push!(dest, j)
		end
	end
	dest
end

findwithin(R::AbstractVector{Vector3D}, Rref::AbstractVector{<:Real},
		dmax::Real) = findwithin!(Int[], R, Rref, nothing, dmax)

findwithin(R::AbstractVector{Vector3D}, Rref::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findwithin!(Int[], R, Rref, cell, dmax)

findwithin(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::AbstractVector{<:Real}, Krefw::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		findwithin!(Int[], Rw, Kw, Rrefw, Krefw, cell, dmax)

findwithin!(dest::AbstractVector{<:Integer}, R::AbstractVector{Vector3D},
		Rref::AbstractVector{<:Real}, dmax::Real) =
		findwithin!(dest, R, Rref, nothing, dmax)

findwithin!(dest::AbstractVector{<:Integer}, R::AbstractVector{Vector3D},
		Rref::AbstractVector{<:Real}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) =
		findwithin!(dest, pbcpos(R, cell)..., pbcpos(Rref, cell)..., cell, dmax)

findwithin!(dest::AbstractVector{<:Integer}, Rw::AbstractVector{Vector3D},
		Kw::AbstractVector{Vector3D}, Rrefw::AbstractVector{<:Real},
		Krefw::AbstractVector{<:Real}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) =
		findwithin!(dest, Rw, Kw, Vector3D(Rrefw), Vector3D(Krefw), cell, dmax)

function findwithin!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real)
	dmax2 = dmax^2
	empty!(dest)
	for i in findnear(Rw, Kw, Rrefw, Krefw, cell, dmax)
		if sqmindist(Rw[i], Kw[i], Rrefw, Krefw, cell) <= dmax2
			push!(dest, i)
		end
	end
	dest
end

struct NeighborDistance
	i::Int
	d2::Real
end

function Base.iterate(nd::NeighborDistance, i::Integer = 1)
	if i > 2
		nothing
	else
		nd[i], i + 1
	end
end

Base.length(::NeighborDistance) = 2

function Base.getindex(nd::NeighborDistance, i::Integer)
	if i == 1
		nd.i
	elseif i == 2
		nd.d2
	else
		throw(BoundsError(nd, i))
	end
end

Base.firstindex(::NeighborDistance) = 1

Base.lastindex(nd::NeighborDistance) = length(nd)

neighborlist(R::AbstractVector{<:Vector3D}, i::Integer, dmax::Real) =
		neighborlist!(NeighborDistance[], R, i, dmax)

neighborlist(R::AbstractVector{<:Vector3D}, i::Integer,
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		neighborlist!(NeighborDistance[], R, i, cell, dmax)

neighborlist(Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		neighborlist!(NeighborDistance[], Rw, Kw, i, cell, dmax)

neighborlist!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		i::Integer, dmax::Real) = neighborlist!(dest, R, i, nothing, dmax)

neighborlist!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		neighborlist!(dest, pbcpos(R, cell)..., i, cell, dmax)

function neighborlist!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real)
	reachable = findnear(Rw, Kw, i, cell, dmax)
	neighborlist!(dest, Rw, Kw, i, cell, dmax, reachable)
end

neighborlist(R::AbstractVector{<:Vector3D}, Rref::AbstractVector{<:Real},
		dmax::Real) = neighborlist!(NeighborDistance[], R, Rref, dmax)

neighborlist(R::AbstractVector{<:Vector3D}, Rref::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		neighborlist!(NeighborDistance[], R, Rref, cell, dmax)

neighborlist(Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		Rrefw::AbstractVector{<:Real}, Krefw::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real) =
		neighborlist!(NeighborDistance[], Rw, Kw, Rrefw, Krefw, cell, dmax)

neighborlist!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		Rref::AbstractVector{<:Real}, dmax::Real) =
		neighborlist!(dest, R, Rref, nothing, dmax)

neighborlist!(dest::AbstractVector{<:Integer}, R::AbstractVector{<:Vector3D},
		Rref::AbstractVector{<:Real}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) = neighborlist!(dest, pbcpos(R, cell)...,
		pbcpos(Rref, cell)..., cell, dmax)

neighborlist!(dest::AbstractVector{<:Integer}, Rw::AbstractVector{<:Vector3D},
		Kw::AbstractVector{<:Vector3D}, Rrefw::AbstractVector{<:Real},
		Krefw::AbstractVector{<:Real}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real) = neighborlist!(dest, Rw, Kw, Vector3D(Rrefw),
		Vector3D(Krefw), cell, dmax)

function neighborlist!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{<:Vector3D}, Kw::AbstractVector{<:Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real)
	reachable = findnear(Rw, Kw, Rrefw, Krefw, cell, dmax)
	neighborlist!(dest, Rw, Kw, Rrefw, Krefw, cell, dmax, reachable)
end

function neighborlist!(dest::AbstractVector{NeighborDistance},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		reachable::AbstractVector{<:Integer})
	dmax2 = dmax^2
	Rwi, Kwi = Rw[i], Kw[i]
	empty!(dest)
	for j in reachable
		if j != i
			d2 = sqmindist(Rw[j], Kw[j], Rwi, Kwi, cell)
			if d2 <= dmax2
				push!(dest, NeighborDistance(j, d2))
			end
		end
	end
	sort!(dest, by = neighbor -> neighbor.d)
end

function neighborlist!(dest::AbstractVector{NeighborDistance},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, reachable::AbstractVector{<:Integer})
	dmax2 = dmax^2
	empty!(dest)
	for i in reachable
		d2 = sqmindist(Rw[i], Kw[i], Rrefw, Krefw, cell)
		if d2 <= dmax2
			push!(dest, NeighborDistance(i, d2))
		end
	end
	sort!(dest, by = neighbor -> neighbor.d)
end

neighborlist(nbrlist::AbstractVector{NeighborDistance}, dmax::Real) =
		neighborlist!(NeighborDistance[], nbrlist, dmax)

function neighborlist!(dest::AbstractVector{NeighborDistance},
		nbrlist::AbstractVector{NeighborDistance}, dmax::Real)
	empty!(dest)
	dmax2 = dmax^2
	i = length(nbrlist)
	while i > 0
		if nbrlist[i].d2 > dmax2
			i -= 1
		else
			break
		end
	end
	sizehint!(dest, length(1:i))
	for j = 1:i
		push!(dest, nbrlist[j])
	end
	dest
end

struct NeighborTableBuffer
	grids::Dict{Int,NeighborGridBuffer}
	lists::Vector{Vector{NeighborDistance}}
	dmaxima::Vector{Float64}
	dest::Vector{Int}
end

NeighborTableBuffer() = NeighborTableBuffer(Dict{Int,NeighborGridBuffer}(),
		Vector{Vector{NeighborDistance}}(), Float64[], Int[])

function prepare!(buffer::NeighborTableBuffer, n::Integer)
	lists = buffer.lists
	dmaxima = buffer.dmaxima
	m = length(dmaxima)
	for list in lists
		empty!(list)
	end
	if n > m
		resize!(lists, n)
		for i = (m+1):n
			lists[i] = NeighborDistance[]
		end
	end
	resize!(dmaxima, n)
	fill!(dmaxima, -1.0)
	buffer
end

struct NeighborTable{T<:NeighborGrid}
	grids::Dict{Int,T}
	gridbuffers::Dict{Int,NeighborGridBuffer}
	lists::Vector{Vector{NeighborDistance}}
	dmaxima::Vector{Float64}
	dest::Vector{Int}

	NeighborTable{T}(grids::Dict{Int,T},
			gridbuffers::Dict{Int,NeighborGridBuffer},
			lists::Vector{Vector{NeighborDistance}}, dmaxima::Vector{Float64},
			dest::Vector{Int}) where {T<:NeighborGrid} =
			new(grids, gridbuffers, lists, dmaxima, dest)
end

function NeighborTable{T}(n::Integer, buffer::NeighborTableBuffer) where
		{T<:NeighborGrid}
	prepare!(buffer, n)
	NeighborTable{T}(Dict{Int,T}(), buffer.grids, buffer.lists, buffer.dmaxima,
			buffer.dest)
end

NeighborTable(n::Integer, cell::TriclinicPBC, buffer::NeighborTableBuffer) =
		NeighborTable{PeriodicNeighborGrid}(n, buffer)

NeighborTable(n::Integer, ::Nothing, buffer::NeighborTableBuffer) =
		NeighborTable{NonperiodicNeighborGrid}(n, buffer)

NeighborTable(n::Integer, cell::Union{TriclinicPBC,Nothing}) =
		NeighborTable(n, cell, NeighborTableBuffer())

Base.length(nbrtable::NeighborTable) = length(nbrtable.dmaxima)

function getnbrgrid(nbrtable::NeighborTable, Rw::AbstractVector{Vector3D},
		Kw::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		dmax::Real)
	dgrid = ceil(Int, dmax)
	get!(nbrtable.grids, dgrid) do
		buffer = get!(nbrtable.gridbuffers, dgrid) do
			NeighborGridBuffer()
		end
		NeighborGrid(Rw, Kw, cell, dmax, buffer)
	end
end

function getnbrlist(nbrtable::NeighborTable,  i::Integer,
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, dmax::Real)
	nbrlist = nbrtable.lists[i]
	if nbrtable.dmaxima[i] < dmax
		nbrgrid = getnbrgrid(nbrtable, Rw, Kw, cell, dmax)
		reachable = findnear!(nbrtable.dest, nbrgrid, i)
		neighborlist!(nbrlist, Rw, Kw, i, cell, dmax, reachable)
	end
	nbrlist
end

findnear(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D}, i::Integer,
		cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable) =
		findnear!(Int[], Rw, Kw, i, cell, dmax, nbrtable)

function findnear!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D}, i::Integer,
		cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable)
	if nbrtable.dmaxima[i] >= dmax
		dmax2 = dmax^2
		empty!(dest)
		push!(dest, i)
		for (j, d2) in nbrtable.lists[i]
			if d2 <= dmax2
				push!(dest, j)
			else
				break
			end
		end
		dest
	else
		nbrgrid = getnbrgrid(nbrtable, Rw, Kw, cell, dmax)
		findnear!(dest, nbrgrid, i)
	end
end

findnear(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable) =
		findnear!(Int[], Rw, Kw, Rrefw, Krefw, cell, dmax, nbrtable)

function findnear!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable)
	nbrgrid = getnbrgrid(nbrtable, Rw, Kw, cell, dmax)
	findnear!(dest, nbrgrid, Krefw)
end

findwithin(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable) =
		findwithin!(Int[], Rw, Kw, i, cell, dmax, nbrtable)

function findwithin!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable)
	dmax2 = dmax^2
	empty!(dest)
	push!(dest, i)
	for (j, d2) in getnbrlist(nbrtable, i, Rw, Kw, cell, dmax)
		if d2 <= dmax2
			push!(dest, i)
		else
			break
		end
	end
	dest
end

findwithin(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable) =
		findwithin!(Int[], Rw, Kw, Rrefw, Krefw, cell, dmax, nbrtable)

function findwithin!(dest::AbstractVector{<:Integer},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable)
	dmax2 = dmax^2
	empty!(dest)
	reachable =
			findnear!(nbrtable.dest, Rw, Kw, Rrefw, Krefw, cell, dmax, nbrtable)
	for i in reachable
		if sqmindist(Rw[i], Kw[i], Rrefw, Krefw, cell) <= dmax2
			push!(dest, i)
		end
	end
	dest
end

neighborlist(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable) =
		neighborlist!(NeighborDistance[], Rw, Kw, i, cell, dmax, nbrtable)

function neighborlist!(dest::AbstractVector{NeighborDistance},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		i::Integer, cell::Union{TriclinicPBC,Nothing}, dmax::Real,
		nbrtable::NeighborTable)
	nbrlist = getnbrlist(nbrtable, i, Rw, Kw, cell, dmax)
	neighborlist!(dest, nbrlist, dmax)
end

neighborlist(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable) = neighborlist!(NeighborDistance[],
		Rw, Kw, Rrefw, Krefw, cell, dmax, nbrtable)

function neighborlist!(dest::AbstractVector{NeighborDistance},
		Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		Rrefw::Vector3D, Krefw::Vector3D, cell::Union{TriclinicPBC,Nothing},
		dmax::Real, nbrtable::NeighborTable)
	reachable =
			findnear!(nbrtable.dest, Rw, Kw, Rrefw, Krefw, cell, dmax, nbrtable)
	neighborlist!(dest, Rw, Kw, Rrefw, Krefw, cell, dmax, reachable)
end

end # module
