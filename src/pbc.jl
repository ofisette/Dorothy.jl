# Periodic boundary conditions

module PBC

using LinearAlgebra
using ..Dorothy.Geometry
using ..Dorothy.Utils
using StaticArrays

export
		TriclinicPBC, RhombododecahedralPBC, OrthorhombicPBC, Kspace, kspace,
		PBCGeometry, pbcgeometry, pbcmatrix, isrhombododecahedral, ishexagonal,
		ismonoclinic, isorthorhombic, istetragonal, iscubic, TriclinicCell,
		RhombododecahedralCell, OrthorhombicCell, pbccell, wrappos, wrappos!,
		pbcpos, pbcpos!, eachimage, nearestpos, nearestpos!, mindist, sqmindist,
		UnwrapStrategy, UnwrapByExtent, UnwrapByGap, unwrap!, PositionGrid,
		posgrid, posgrid!

abstract type TriclinicPBC end

abstract type RhombododecahedralPBC <: TriclinicPBC end

abstract type OrthorhombicPBC <: TriclinicPBC end

struct Kspace end

const kspace = Kspace()

struct PBCGeometry
	sides::Vector3D
	angles::Vector3D
end

function Base.iterate(geometry::PBCGeometry, state::Integer = 1)
	if state == 1
		geometry.sides, 2
	elseif state == 2
		geometry.angles, 3
	else
		nothing
	end
end

Base.eltype(::Type{PBCGeometry}) = Vector3D

Base.length(::PBCGeometry) = 2

function Base.getindex(geometry::PBCGeometry, i::Integer)
	if i == 1
		geometry.sides
	elseif i == 2
		geometry.angles
	else
		throw(BoundsError(geometry, i))
	end
end

pbcgeometry(cell::TriclinicPBC) = pbcgeometry(pbcmatrix(cell))

pbcgeometry(M::AbstractMatrix{<:Real}) = pbcgeometry(Basis3D(M))

function pbcgeometry(M::Basis3D)
	A, B, C = [M[:,i] for i = 1:3]
	a, b, c = [norm(V) for V in (A, B, C)]
	α, β, γ = acos(B⋅C / (b*c)), acos(A⋅C / (a*c)), acos(A⋅B / (a*b))
	PBCGeometry([a,b,c], [α,β,γ])
end

pbcgeometry(sides::AbstractVector{<:Real}, angles::AbstractVector{<:Real}) =
		PBCGeometry(sides, angles)

pbcgeometry(geometry::PBCGeometry) = geometry

pbcmatrix(M::AbstractMatrix{<:Real}) = Basis3D(M)

pbcmatrix(M::Basis3D) = M

pbcmatrix(sides::AbstractVector{<:Real}, angles::AbstractVector{<:Real}) =
		pbcmatrix(pbcgeometry(sides, angles))

function pbcmatrix(((a,b,c),(α,β,γ))::PBCGeometry)
	Ax, Ay, Az = a, 0.0, 0.0 # A lies along the positive X axis
	Bx, By, Bz = b*cos(α), b*sin(α), 0.0 # B lies in the X-Y plane
	Cx = c*cos(β)
	Cy = (b*c*cos(γ) - Bx*Cx) / By
	Cz = √(c^2 - Cx^2 - Cy^2)
	[[Ax,Ay,Az] [Bx,By,Bz] [Cx,Cy,Cz]]
end

Base.:(==)(cell1::TriclinicPBC, cell2::TriclinicPBC) =
		pbcmatrix(cell1) == pbcmatrix(cell2)

Base.:(*)(cell::TriclinicPBC, R::AbstractVector{<:Real}) = cell.T * R

Base.broadcastable(cell::TriclinicPBC) = Ref(cell)

Base.inv(cell::TriclinicPBC) = cell.inv

Geometry.dims(cell::TriclinicPBC) = Vector3D(sum(pbcmatrix(cell); dims=2))

Geometry.center(cell::TriclinicPBC) = dims(cell) ./ 2

Geometry.extent(cell::TriclinicPBC) =
		BoundingBox(Vector3D(0.0, 0.0, 0.0), dims(cell))

Geometry.isinside(R::AbstractVector{<:Real}, cell::TriclinicPBC) =
		isinside(Vector3D(R), cell)

Geometry.isinside(R::Vector3D, cell::TriclinicPBC) =
		isinside(inv(cell) * R, kspace)

Geometry.isinside(R::Vector3D, cell::OrthorhombicPBC) =
		(0.0 <= R.x < dims(cell).x) && (0.0 <= R.y < dims(cell).y) &&
		(0.0 <= R.z < dims(cell).z)

Geometry.isinside(Rk::Vector3D, ::Kspace) =
		(0.0 <= Rk.x < 1.0) && (0.0 <= Rk.y < 1.0) && (0.0 <= Rk.z < 1.0)

function Geometry.volume(((a,b,c),(α,β,γ))::PBCGeometry)
	cosα, cosβ, cosγ = cos(α), cos(β), cos(γ)
	a*b*c * √(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
end

function Geometry.volume(cell::TriclinicPBC)
	M = pbcmatrix(cell)
	A, B, C = [M[1:3,i] for i = 1:3]
	A⋅(B×C)
end

isrhombododecahedral(((a,b,c),(α,β,γ))::PBCGeometry) =
		(a ≈ b ≈ c) && (α ≈ τ/4) && (β ≈ γ ≈ τ/6)

isrhombododecahedral(cell::TriclinicPBC) =
		isrhombododecahedral(pbcgeometry(cell))

isrhombododecahedral(cell::RhombododecahedralPBC) = true

isrhombododecahedral(cell::OrthorhombicPBC) = false

ishexagonal(((a,b,c),(α,β,γ))::PBCGeometry) =
		(a ≈ b) && (α ≈ β ≈ τ/4) && (γ ≈ τ/3)

ishexagonal(cell::TriclinicPBC) = ishexagonal(pbcgeometry(cell))

ishexagonal(cell::RhombododecahedralPBC) = false

ishexagonal(cell::OrthorhombicPBC) = false

ismonoclinic(((a,b,c),(α,β,γ))::PBCGeometry) = (α ≈ τ/4)

ismonoclinic(cell::TriclinicPBC) = ismonoclinic(pbcgeometry(cell))

ismonoclinic(cell::RhombododecahedralPBC) = false

ismonoclinic(cell::OrthorhombicPBC) = true

isorthorhombic(((a,b,c),(α,β,γ))::PBCGeometry) = (α ≈ β ≈ γ ≈ τ/4)

isorthorhombic(cell::TriclinicPBC) = isorthorhombic(pbcgeometry(cell))

isorthorhombic(cell::RhombododecahedralPBC) = false

isorthorhombic(cell::OrthorhombicPBC) = true

istetragonal(((a,b,c),(α,β,γ))::PBCGeometry) = (a ≈ b) && (α ≈ β ≈ γ ≈ τ/4)

istetragonal(cell::TriclinicPBC) = istetragonal(pbcgeometry(cell))

istetragonal(cell::RhombododecahedralPBC) = false

istetragonal(cell::OrthorhombicPBC) = (dims(cell).x ≈ dims(cell).y)

iscubic(((a,b,c),(α,β,γ))::PBCGeometry) = (a ≈ b ≈ c) && (α ≈ β ≈ γ ≈ τ/4)

iscubic(cell::TriclinicPBC) = iscubic(pbcgeometry(cell))

iscubic(cell::RhombododecahedralPBC) = false

iscubic(cell::OrthorhombicPBC) = (dims(cell).x ≈ dims(cell).y ≈ dims(cell).z)

struct TriclinicCell <: TriclinicPBC
	T::LinearTransformation
	inv::LinearTransformation
	r2::Float64

	function TriclinicCell(M::AbstractMatrix{<:Real})
		r2 = (minimum(pbcgeometry(M).sides) / 2)^2
		new(LinearTransformation(M), LinearTransformation(inv(M)), r2)
	end
end

TriclinicCell(sides::AbstractVector{<:Real}, angles::AbstractVector{<:Real}) =
		TriclinicCell(pbcmatrix(sides, angles))

TriclinicCell(geometry::PBCGeometry) = TriclinicCell(pbcmatrix(geometry))

TriclinicCell(cell::TriclinicPBC) = TriclinicCell(pbcmatrix(cell))

(cell::TriclinicCell)(x) = transform(cell.T, x)

pbcmatrix(cell::TriclinicCell) = cell.T.M

struct RhombododecahedralCell <: RhombododecahedralPBC
	T::LinearTransformation
	inv::LinearTransformation
	r2::Float64

	function RhombododecahedralCell(a::Real)
		M = pbcmatrix([a,a,a], [τ/4, τ/6, τ/6])
		new(LinearTransformation(M), LinearTransformation(inv(M)), (a/2)^2)
	end
end

RhombododecahedralCell(cell::RhombododecahedralCell) = cell

pbcmatrix(cell::RhombododecahedralCell) = cell.T.M

pbcgeometry(cell::RhombododecahedralCell) =
		(dims(cell), Vector3D(τ/4, τ/6, τ/6))

(cell::RhombododecahedralCell)(x) = transform(cell.T, x)

struct OrthorhombicCell <: OrthorhombicPBC
	T::Scaling
	inv::Scaling
	R::Vector3D

	OrthorhombicCell(V::AbstractVector{<:Real}) =
			new(Scaling(V), Scaling(1.0 ./ V), V ./ 2)
end

OrthorhombicCell(xy::Real, z::Real) = OrthorhombicCell([xy, xy, z])

OrthorhombicCell(xyz::Real) = OrthorhombicCell([xyz, xyz, xyz])

OrthorhombicCell(cell::OrthorhombicCell) = cell

pbcmatrix(cell::OrthorhombicCell) = Diagonal(cell.T.V)

pbcgeometry(cell::OrthorhombicCell) = (cell.T.V, Vector3D(τ/4, τ/4, τ/4))

Geometry.dims(cell::OrthorhombicCell) = cell.T.V

Geometry.center(cell::OrthorhombicCell) = cell.R

(cell::OrthorhombicCell)(x) = transform(cell.T, x)

pbccell(cell::TriclinicPBC) = pbccell(pbcmatrix(cell))

pbccell(cell::RhombododecahedralCell) = cell

pbccell(cell::OrthorhombicCell) = cell

pbccell(M::AbstractMatrix{<:Real}) = pbccell(Basis3D(M))

function pbccell(M::Basis3D)
	geometry = pbcgeometry(M)
	if isrhombododecahedral(geometry)
		RhombododecahedralCell(geometry.sides[1])
	elseif isorthorhombic(geometry)
		OrthorhombicCell(geometry.sides)
	else
		TriclinicCell(M)
	end
end

pbccell(sides::AbstractVector{<:Real}, angles::AbstractVector{<:Real}) =
		pbccell(PBCGeometry(sides, angles))

function pbccell(geometry::PBCGeometry)
	if isrhombododecahedral(geometry)
		RhombododecahedralCell(geometry.sides[1])
	elseif isorthorhombic(geometry)
		OrthorhombicCell(geometry.sides)
	else
		TriclinicCell(pbcmatrix(geometry))
	end
end

pbccell(::Nothing) = nothing

wrappos(R::AbstractVector{<:Real}, cell::TriclinicPBC) =
		wrappos(Vector3D(R), cell)

wrappos(R::Vector3D, cell::TriclinicPBC) = cell * wrappos(inv(cell) * R, kspace)

wrappos(Rk::Vector3D, ::Kspace) = (Rk - floor.(Rk))

wrappos!(R::AbstractVector{Vector3D}, cell::TriclinicPBC) =
		wrappos!(cog, R, cell)

function wrappos!(f::Function, R::AbstractVector{Vector3D}, cell::TriclinicPBC)
	@boundscheck length(R) > 0 || error("cannot wrap zero positions")
	C = f(R)
	T = wrappos(C, cell) - C
	if ! iszero(T)
		translation(T)(R)
	end
	R
end

pbcpos(R::AbstractVector{<:Real}, cell::TriclinicPBC) =
		pbcpos(Vector3D(R), cell)

function pbcpos(R::Vector3D, cell::TriclinicPBC)
	Rkw = wrappos(inv(cell) * R, kspace)
	Rw = cell * Rkw
	Rw, Rkw
end

pbcpos(R::Vector3D, ::Nothing) = R, R

pbcpos(R::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing}) =
		pbcpos!(similar(R), similar(R), R, cell)

function pbcpos!(Rw::AbstractVector{Vector3D}, Rkw::AbstractVector{Vector3D},
		R::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing})
	@boundscheck length(Rw) == length(Rkw) == length(R) ||
			error("size mismatch between position arrays")
	for i in eachindex(R)
		Rw[i], Rkw[i] = pbcpos(R[i], cell)
	end
	Rw, Rkw
end

struct PBCImageIterator{T<:TriclinicPBC}
	Rkw::Vector3D
	cell::T

	PBCImageIterator{T}(Rkw::Vector3D, cell::T) where {T<:TriclinicPBC} =
			new(Rkw, cell)
end

PBCImageIterator(Rkw::Vector3D, cell::T) where {T<:TriclinicPBC} =
		PBCImageIterator{T}(Rkw, cell)

Base.eltype(::Type{PBCImageIterator}) = Vector3D

Base.length(::PBCImageIterator) = 26

const pbcfirstshell = @SVector [
	Vector3D(-1.0,-1.0,-1.0), Vector3D(-1.0,-1.0,0.0), Vector3D(-1.0,-1.0,1.0),
	Vector3D(-1.0,0.0,-1.0),  Vector3D(-1.0,0.0,0.0),  Vector3D(-1.0,0.0,1.0),
	Vector3D(-1.0,1.0,-1.0),  Vector3D(-1.0,1.0,0.0),  Vector3D(-1.0,1.0,1.0),
	Vector3D(0.0,-1.0,-1.0),  Vector3D(0.0,-1.0,0.0),  Vector3D(0.0,-1.0,1.0),
	Vector3D(0.0,0.0,-1.0),                            Vector3D(0.0,0.0,1.0),
	Vector3D(0.0,1.0,-1.0),   Vector3D(0.0,1.0,0.0),   Vector3D(0.0,1.0,1.0),
	Vector3D(1.0,-1.0,-1.0),  Vector3D(1.0,-1.0,0.0),  Vector3D(1.0,-1.0,1.0),
	Vector3D(1.0,0.0,-1.0),   Vector3D(1.0,0.0,0.0),   Vector3D(1.0,0.0,1.0),
	Vector3D(1.0,1.0,-1.0),   Vector3D(1.0,1.0,0.0),   Vector3D(1.0,1.0,1.0)]

function Base.iterate(iter::PBCImageIterator, n = 1)
	if n > 26
		nothing
	else
		iter.cell * (iter.Rkw + pbcfirstshell[n]), n + 1
	end
end

eachimage(R::Vector3D, cell::TriclinicPBC) =
		PBCImageIterator(wrappos(inv(cell) * R, kspace), cell)

eachimage(Rw::Vector3D, Rkw::Vector3D, cell::TriclinicPBC) =
		PBCImageIterator(Rkw, cell)

nearestpos(R::AbstractVector{<:Real}, O::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		nearestpos(Vector3D(R), Vector3D(O), cell)

nearestpos(R::Vector3D, O::Vector3D, cell::Union{TriclinicPBC,Nothing}) =
		nearestpos(pbcpos(R, cell)..., O, cell)

function nearestpos(Rw::Vector3D, Rkw::Vector3D, O::Vector3D,
		cell::TriclinicPBC)
	Rwmin = Rw
	d2min = sqdist(Rw, O)
	if d2min < cell.r2
		return Rwmin
	end
	for Rwi in eachimage(Rw, Rkw, cell)
		d2 = sqdist(Rwi, O)
		if d2 < cell.r2
			return Rwi
		elseif d2 < d2min
			d2min = d2
			Rwmin = Rwi
		end
	end
	Rwmin
end

function nearestpos(Rw::Vector3D, Rkw::Vector3D, O::Vector3D,
		cell::OrthorhombicPBC)
	(dx,dy,dz) = O - Rw
	(x,y,z) = Rw
	if dx > center(cell).x
		x += dims(cell).x
	elseif dx < -center(cell).x
		x -= dims(cell).x
	end
	if dy > center(cell).y
		y += dims(cell).y
	elseif dy < -center(cell).y
		y -= dims(cell).y
	end
	if dz > center(cell).z
		z += dims(cell).z
	elseif dz < -center(cell).z
		z -= dims(cell).z
	end
	Vector3D(x,y,z)
end

nearestpos(Rw::Vector3D, Rkw::Vector3D, O::Vector3D, ::Nothing) = Rw

nearestpos!(R::AbstractVector{Vector3D}, O::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) = nearestpos!(R, Vector3D(O), cell)

nearestpos!(f, R::AbstractVector{Vector3D}, O::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		nearestpos!(f, R, Vector3D(O), cell)

nearestpos!(R::AbstractVector{Vector3D}, O::Vector3D,
		cell::Union{TriclinicPBC,Nothing}) = nearestpos!(cog, R, O, cell)

function nearestpos!(f, R::AbstractVector{Vector3D}, O::Vector3D,
		cell::Union{TriclinicPBC,Nothing})
	@boundscheck length(R) > 0 || error("cannot translate zero positions")
	C::Vector3D = f(R)
	T = nearestpos(C, O, cell) - C
	if ! iszero(T)
		translation(T)(R)
	end
	R
end

mindist(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) = sqrt(sqmindist(R1, R2, cell))

mindist(R1w::AbstractVector{<:Real}, R1kw::AbstractVector{<:Real},
		R2w::AbstractVector{<:Real}, R2kw::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		sqrt(sqmindist(R1w, R1kw, R2w, R2kw, cell))

sqmindist(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		sqmindist(Vector3D(R1), Vector3D(R2), cell)

sqmindist(R1w::AbstractVector{<:Real}, R1kw::AbstractVector{<:Real},
		R2w::AbstractVector{<:Real}, R2kw::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) = sqmindist(Vector3D(R1w),
		Vector3D(R1kw), Vector3D(R2w), Vector3D(R2kw), cell)

sqmindist(R1::Vector3D, R2::Vector3D, cell::Union{TriclinicPBC,Nothing}) =
		sqmindist(pbcpos(R1)..., pbcpos(R2)..., cell)

function sqmindist(R1w::Vector3D, R1kw::Vector3D, R2w::Vector3D, R2kw::Vector3D,
		cell::TriclinicPBC)
	d2min = sqdist(R1w, R2w)
	if d2min < cell.r2
		return d2min
	end
	for R1wi in eachimage(R1w, R1kw, cell)
		d2 = sqdist(R1wi, R2w)
		if d2 < cell.r2
			return d2
		elseif d2 < d2min
			d2min = d2
		end
	end
	d2min
end

function sqmindist(R1w::Vector3D, R1kw::Vector3D, R2w::Vector3D, R2kw::Vector3D,
		cell::OrthorhombicPBC)
	x, y, z = abs.(R2w - R1w)
	if x > center(cell).x
		x = dims(cell).x - x
	end
	if y > center(cell).y
		y = dims(cell).y - y
	end
	if z > center(cell).z
		z = dims(cell).z - z
	end
	x^2 + y^2 + z^2
end

sqmindist(R1w::Vector3D, R1kw::Vector3D, R2w::Vector3D, R2kw::Vector3D,
		::Nothing) = sqdist(R1w, R2w)

abstract type UnwrapStrategy end

struct UnwrapByExtent <: UnwrapStrategy end

function unwrap!(R::AbstractVector{Vector3D}, cell::TriclinicPBC,
		strategy::UnwrapByExtent)
	@boundscheck length(R) > 1 || error("cannot unwrap fewer than 2 positions")
	inv(cell)(R)
	unwrap!(R, kspace, strategy)
	cell(R)
end

function unwrap!(Rk::AbstractVector{Vector3D}, ::Kspace,
		strategy::UnwrapByExtent)
	@boundscheck length(Rk) > 1 || error("cannot unwrap fewer than 2 positions")
	Rkdims = dims(extent(Rk))
	if Rkdims.x > 0.5
		for i in eachindex(Rk)
			x, y, z = Rk[i]
			x -= 0.5
			x -= floor(x)
			x += 0.5
			Rk[i] = Vector3D(x,y,z)
		end
	end
	if Rkdims.y > 0.5
		for i in eachindex(Rk)
			x, y, z = Rk[i]
			y -= 0.5
			y -= floor(y)
			y += 0.5
			Rk[i] = Vector3D(x,y,z)
		end
	end
	if Rkdims.z > 0.5
		for i in eachindex(Rk)
			x, y, z = Rk[i]
			z -= 0.5
			z -= floor(z)
			z += 0.5
			Rk[i] = Vector3D(x,y,z)
		end
	end
	Rk
end

struct UnwrapByGap <: UnwrapStrategy
	D::Vector3D
	buffer::Vector{Vector3D}

	function UnwrapByGap(cell::TriclinicPBC, d::Real = 2.5)
		(a,b,c), angles = pbcgeometry(cell)
		new(Vector3D(d,d,d) ./ Vector3D(a,b,c), Vector3D[])
	end
end

unwrap!(R::AbstractVector{Vector3D}, cell::TriclinicPBC) =
		unwrap!(R, cell, UnwrapByGap(cell))

function unwrap!(R::AbstractVector{Vector3D}, cell::TriclinicPBC,
		strategy::UnwrapByGap)
	@boundscheck length(R) > 1 || error("cannot unwrap fewer than 2 positions")
	inv(cell)(R)
	unwrap!(R, kspace, strategy)
	cell(R)
end

function unwrap!(Rk::AbstractVector{Vector3D}, ::Kspace,
		strategy::UnwrapByGap)
	@boundscheck length(Rk) > 1 || error("cannot unwrap fewer than 2 positions")
	Rks = resize!(strategy.buffer, length(Rk))
	Rks .= Rk
	sort!(Rks, by = R -> R.x)
	for i = 2:length(Rks)
		if Rks[i].x - Rks[i-1].x > strategy.D.x
			d = Rks[i].x
			for j in eachindex(Rk)
				x, y, z = Rk[j]
				x -= d
				x -= floor(x)
				x += d
				Rk[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	sort!(Rks, by = R -> R.y)
	for i = 2:length(Rks)
		if Rks[i].y - Rks[i-1].y > strategy.D.y
			d = Rks[i].y
			for j in eachindex(Rk)
				x, y, z = Rk[j]
				y -= d
				y -= floor(y)
				y += d
				Rk[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	sort!(Rks, by = R -> R.z)
	for i = 2:length(Rks)
		if Rks[i].z - Rks[i-1].z > strategy.D.z
			d = Rks[i].z
			for j in eachindex(Rk)
				x, y, z = Rk[j]
				z -= d
				z -= floor(z)
				z += d
				Rk[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	Rk
end

abstract type PositionGrid end

struct NonperiodicPositionGrid <: PositionGrid
	g3::NonperiodicGrid3D{Int}
	I::Vector{Tuple{Int,Int,Int}}

	NonperiodicPositionGrid(g3::NonperiodicGrid3D{Int},
			I::AbstractVector{<:Tuple{Integer,Integer,Integer}}) = new(g3, I)
end

function NonperiodicPositionGrid()
	g3 = NonperiodicGrid3D{Int}(extent([0.0,0.0,0.0]), [1.0,1.0,1.0])
	I = Vector{Tuple{Int,Int,Int}}(undef, 0)
	NonperiodicPositionGrid(g3, I)
end

NonperiodicPositionGrid(R::AbstractVector{Vector3D}, d::Real) =
		NonperiodicPositionGrid(NonperiodicPositionGrid(), R, d)

function NonperiodicPositionGrid(pg::NonperiodicPositionGrid,
		R::AbstractVector{Vector3D}, d::Real)
	g3 = resize!(pg.g3, extent(R), [d,d,d])
	I = resize!(pg.I, length(R))
	for i in eachindex(R)
		(x,y,z) = findcell(g3, R[i])
		push!(g3, (x,y,z), i)
		I[i] = (x,y,z)
	end
	NonperiodicPositionGrid(g3, I)
end

struct PeriodicPositionGrid{T<:TriclinicPBC} <: PositionGrid
	g3::PeriodicGrid3D{Int}
	I::Vector{Tuple{Int,Int,Int}}
	cell::T

	PeriodicPositionGrid{T}(g3::PeriodicGrid3D{Int},
			I::AbstractVector{<:Tuple{Integer,Integer,Integer}}, cell::T) where
			{T<:TriclinicPBC} = new(g3, I, cell)
end

PeriodicPositionGrid(g3::PeriodicGrid3D{Int},
		I::AbstractVector{<:Tuple{Integer,Integer,Integer}}, cell::T) where
		{T<:TriclinicPBC} = PeriodicPositionGrid{T}(g3, I, cell)

function PeriodicPositionGrid(cell::TriclinicPBC)
	g3 = PeriodicGrid3D{Int}(extent([0.0,0.0,0.0]), [1.0,1.0,1.0])
	I = Vector{Tuple{Int,Int,Int}}(undef, 0)
	PeriodicPositionGrid(g3, I, cell)
end

PeriodicPositionGrid(Rkw::AbstractVector{Vector3D}, cell::TriclinicPBC,
		d::Real) =
		PeriodicPositionGrid(PeriodicPositionGrid(cell), Rkw, cell, d)

function PeriodicPositionGrid(pg::PeriodicPositionGrid,
		Rkw::AbstractVector{Vector3D}, cell::TriclinicPBC, d::Real)
	D = [1.0,1.0,1.0] - inv(cell) * (dims(cell) .- d)
	g3 = resize!(pg.g3, extent([0.0,0.0,0.0], [1.0,1.0,1.0]), D)
	I = resize!(pg.I, length(Rkw))
	for i in eachindex(Rkw)
		(x,y,z) = findcell(g3, Rkw[i])
		push!(g3, (x,y,z), i)
		I[i] = (x,y,z)
	end
	PeriodicPositionGrid(g3, I, cell)
end

Geometry.findnear(pg::PositionGrid, i::Integer) = findnear!(Int[], pg, i)

Geometry.findnear(pg::PositionGrid, R::AbstractVector{<:Real}) =
		findnear!(Int[], pg, R)

Geometry.findnear(pg::PositionGrid, Rw::AbstractVector{<:Real},
		Rkw::AbstractVector{<:Real}) = findnear!(Int[], pg, Rw, Rkw)

Geometry.findnear!(dest::AbstractVector{<:Integer}, pg::PositionGrid,
		i::Integer) = findnear!(dest, pg.g3, pg.I[i])

Geometry.findnear!(dest::AbstractVector{<:Integer}, pg::PositionGrid,
		R::AbstractVector{<:Real}) = findnear!(dest, pg, pbcpos(R, pg.cell)...)

Geometry.findnear!(dest::AbstractVector{<:Integer}, pg::PositionGrid,
		Rw::AbstractVector{<:Real}, Rkw::AbstractVector{<:Real}) =
		findnear!(dest, pg, Vector3D(Rw), Vector3D(Rkw))

Geometry.findnear!(dest::AbstractVector{<:Integer}, pg::PositionGrid,
		Rw::Vector3D, Rkw::Vector3D) = findnear!(dest, pg.g3, Rkw)

posgrid(::Nothing) = NonperiodicPositionGrid()

posgrid(cell::TriclinicPBC) = PeriodicPositionGrid(cell)

posgrid(R::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing},
		d::Real) = posgrid(pbcpos(R, cell)..., cell, d)

posgrid(Rw::AbstractVector{Vector3D}, Rkw::AbstractVector{Vector3D},
		::Nothing, d::Real) = NonperiodicPositionGrid(Rw, d)

posgrid(Rw::AbstractVector{Vector3D}, Rkw::AbstractVector{Vector3D},
		cell::TriclinicPBC, d::Real) = PeriodicPositionGrid(Rkw, cell, d)

posgrid!(pg::PositionGrid, R::AbstractVector{Vector3D},
		cell::Union{TriclinicPBC,Nothing}, d::Real) =
		posgrid(pg, pbcpos(R, cell)..., cell, d)

posgrid!(pg::NonperiodicPositionGrid, Rw::AbstractVector{Vector3D},
		Rkw::AbstractVector{Vector3D}, ::Nothing, d::Real) =
		NonperiodicPositionGrid(pg, Rw, d)

posgrid!(pg::PeriodicPositionGrid, Rw::AbstractVector{Vector3D},
		Rkw::AbstractVector{Vector3D}, cell::TriclinicPBC, d::Real) =
		PeriodicPositionGrid(pg, Rkw, cell, d)

end # module
