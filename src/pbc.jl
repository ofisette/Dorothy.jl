# Periodic boundary conditions

module PBC

using LinearAlgebra
using ..Dorothy.Geometry
using ..Dorothy.Properties
using ..Dorothy.Utils
using StaticArrays

export
		TriclinicPBC, RhombododecahedralPBC, OrthorhombicPBC, Kspace, kspace,

		pbcgeometry, pbcmatrix, pbcvolume, isrhombododecahedral,
		istruncatedoctahedral, ishexagonal, isorthorhombic, istetragonal,
		iscubic,

		TriclinicCell, RhombododecahedralCell, OrthorhombicCell, pbccell,

		wrappos, wrappos!, pbcpos, pbcpos!, eachimage, nearestpos, nearestpos!,
		mindist, sqmindist, UnwrapStrategy, UnwrapByExtent, UnwrapByGap,
		unwrap!,

		ProximityLattice, NonperiodicProximityLattice, PeriodicProximityLattice,
		proxiparams, proxilattice

const RealTriple = Tuple{Real,Real,Real}

abstract type TriclinicPBC end

abstract type RhombododecahedralPBC <: TriclinicPBC end

abstract type OrthorhombicPBC <: TriclinicPBC end

struct Kspace end

const kspace = Kspace()

pbcgeometry(cell::TriclinicPBC) = pbcgeometry(pbcmatrix(cell))

pbcgeometry(M::AbstractMatrix{<:Real}) = pbcgeometry(Basis3D(M))

function pbcgeometry(M::Basis3D)
	A, B, C = [M[:,i] for i = 1:3]
	a, b, c = [norm(V) for V in (A, B, C)]
	α, β, γ = acos(B⋅C / (b*c)), acos(A⋅C / (a*c)), acos(A⋅B / (a*b))
	(a,b,c), (α,β,γ)
end

function pbcmatrix((a,b,c)::RealTriple, (α,β,γ)::RealTriple)
	Ax, Ay, Az = a, 0.0, 0.0 # A lies along the positive X axis
	Bx, By, Bz = b*cos(γ), b*sin(γ), 0.0 # B lies in the X-Y plane
	Cx = c*cos(β)
	Cy = (b*c*cos(α) - Bx*Cx) / By
	Cz = √(c^2 - Cx^2 - Cy^2)
	Basis3D([Ax Bx Cx;
	         Ay By Cy;
			 Az Bz Cz])
end

Base.:(==)(cell1::TriclinicPBC, cell2::TriclinicPBC) =
		pbcmatrix(cell1) ≈ pbcmatrix(cell2)

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

Geometry.isinside(K::Vector3D, ::Kspace) =
		(0.0 <= K.x < 1.0) && (0.0 <= K.y < 1.0) && (0.0 <= K.z < 1.0)

function pbcvolume((a,b,c)::RealTriple, (α,β,γ)::RealTriple)
	cosα, cosβ, cosγ = cos(α), cos(β), cos(γ)
	a*b*c * √(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
end

pbcvolume(cell::TriclinicPBC) = pbcvolume(pbcmatrix(cell))

pbcvolume(M::AbstractMatrix{<:Real}) = pbcvolume(Basis3D(M))

function pbcvolume(M::Basis3D)
	A, B, C = [M[1:3,i] for i = 1:3]
	A⋅(B×C)
end

Geometry.volume(cell::TriclinicPBC) = pbcvolume(cell)

pbcrms(((a1,b1,c1),(α1,β1,γ1))::Tuple{RealTriple,RealTriple},
		((a2,b2,c2),(α2,β2,γ2))::Tuple{RealTriple,RealTriple}) =
		sqrt(abs(a2-a1)^2 + abs(b2-b1)^2 + abs(c2-c1)^2),
		sqrt(abs(α2-α1)^2 + abs(β2-β1)^2 + abs(γ2-γ1)^2)

const defedgetol = 0.01
const defangletol = 0.01

function isequalpbc(((a1,b1,c1),(α1,β1,γ1))::Tuple{RealTriple,RealTriple},
		((a2,b2,c2),(α2,β2,γ2))::Tuple{RealTriple,RealTriple};
		edgetol::Real = defedgetol, angletol::Real = defangletol)
	rmsedges, rmsangles =
			pbcrms(((a1,b1,c1),(α1,β1,γ1)), ((a2,b2,c2),(α2,β2,γ2)))
	rmsedges <= edgetol && rmsangles <= angletol
end

rhombododecahedral((a,b,c)::RealTriple, (α,β,γ)::RealTriple) =
		(a, a, a), (τ/6, τ/6, τ/4)

isrhombododecahedral((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), rhombododecahedral((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

function truncatedoctahedral((a,b,c)::RealTriple,(α,β,γ)::RealTriple)
	θ = arcos(1/3)
	(a, a, a), (θ, τ/2 - θ, θ)
end

istruncatedoctahedral((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), truncatedoctahedral((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

hexagonal((a,b,c)::RealTriple, (α,β,γ)::RealTriple) =
		(a, a, c), (τ/4, τ/4, τ/3)

ishexagonal((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), hexagonal((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

orthorhombic((a,b,c)::RealTriple, (α,β,γ)::RealTriple) =
		(a, b, c), (τ/4, τ/4, τ/4)

isorthorhombic((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), orthorhombic((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

tetragonal((a,b,c)::RealTriple,(α,β,γ)::RealTriple) = (a, a, c), (τ/4, τ/4, τ/4)

istetragonal((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), tetragonal((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

cubic((a,b,c)::RealTriple, (α,β,γ)::RealTriple) = (a, a, a), (τ/4, τ/4, τ/4)

iscubic((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol) =
		isequalpbc(((a,b,c), (α,β,γ)), cubic((a,b,c), (α,β,γ));
		edgetol=edgetol, angletol=angletol)

struct TriclinicCell <: TriclinicPBC
	T::LinearTransformation
	inv::LinearTransformation

	TriclinicCell(M::AbstractMatrix{<:Real}) =
			new(LinearTransformation(M), LinearTransformation(inv(M)))
end

TriclinicCell((a,b,c)::RealTriple,(α,β,γ)::RealTriple) =
		TriclinicCell(pbcmatrix((a,b,c), (α,β,γ)))

TriclinicCell(cell::TriclinicPBC) = TriclinicCell(pbcmatrix(cell))

(cell::TriclinicCell)(x) = transform(cell.T, x)

pbcmatrix(cell::TriclinicCell) = cell.T.M

struct RhombododecahedralCell <: RhombododecahedralPBC
	T::LinearTransformation
	inv::LinearTransformation
	r2::Float64

	function RhombododecahedralCell(a::Real)
		M = pbcmatrix((a,a,a), (τ/6, τ/6, τ/4))
		new(LinearTransformation(M), LinearTransformation(inv(M)), (a/2)^2)
	end
end

pbcmatrix(cell::RhombododecahedralCell) = cell.T.M

function pbcgeometry(cell::RhombododecahedralCell)
	a = 2*sqrt(cell.r2)
	(a,a,a), (τ/6,τ/6,τ/4)
end

(cell::RhombododecahedralCell)(x) = transform(cell.T, x)

struct OrthorhombicCell <: OrthorhombicPBC
	T::Scaling
	inv::Scaling

	OrthorhombicCell(V::AbstractVector{<:Real}) =
			new(Scaling(V), Scaling(1.0 ./ V))
end

OrthorhombicCell(x::Real, y::Real, z::Real) = OrthorhombicCell([x, y, z])

OrthorhombicCell(xy::Real, z::Real) = OrthorhombicCell([xy, xy, z])

OrthorhombicCell(xyz::Real) = OrthorhombicCell([xyz, xyz, xyz])

pbcmatrix(cell::OrthorhombicCell) = Diagonal(cell.T.V)

pbcgeometry(cell::OrthorhombicCell) = (dims(cell)...,), (τ/4,τ/4,τ/4)

Geometry.dims(cell::OrthorhombicCell) = cell.T.V

Geometry.volume(cell::OrthorhombicCell) = prod(dims(cell))

(cell::OrthorhombicCell)(x) = transform(cell.T, x)

pbccell(::Nothing; edgetol::Real = defedgetol,
		angletol::Real = defangletol) = nothing

pbccell(cell::TriclinicPBC; edgetol::Real = defedgetol,
		angletol::Real = defangletol) =
		pbccell(pbcgeometry(cell)...; edgetol=edgetol, angletol=angletol)

pbccell(cell::RhombododecahedralCell; edgetol::Real = defedgetol,
		angletol::Real = defangletol) = cell

pbccell(cell::OrthorhombicCell; edgetol::Real = defedgetol,
		angletol::Real = defangletol) = cell

pbccell(M::AbstractMatrix{<:Real}; edgetol::Real = defedgetol,
		angletol::Real = defangletol) =
		pbccell(pbcgeometry(M)...; edgetol=edgetol, angletol=angletol)

function pbccell((a,b,c)::RealTriple, (α,β,γ)::RealTriple;
		edgetol::Real = defedgetol, angletol::Real = defangletol)
	if isrhombododecahedral((a,b,c), (α,β,γ); edgetol=edgetol,
			angletol=angletol)
		RhombododecahedralCell(a)
	elseif isorthorhombic((a,b,c), (α,β,γ); edgetol=edgetol,
			angletol=angletol)
		OrthorhombicCell(a, b, c)
	else
		TriclinicCell((a,b,c), (α,β,γ))
	end
end

wrappos(R::AbstractVector{<:Real}, cell::TriclinicPBC) =
		wrappos(Vector3D(R), cell)

wrappos(R::Vector3D, cell::TriclinicPBC) = cell * wrappos(inv(cell) * R, kspace)

wrappos(K::Vector3D, ::Kspace) = (K - floor.(K))

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
	Kw = wrappos(inv(cell) * R, kspace)
	Rw = cell * Kw
	Rw, Kw
end

pbcpos(R::Vector3D, ::Nothing) = R, R

pbcpos(R::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing}) =
		pbcpos!(similar(R), similar(R), R, cell)

function pbcpos!(Rw::AbstractVector{Vector3D}, Kw::AbstractVector{Vector3D},
		R::AbstractVector{Vector3D}, cell::Union{TriclinicPBC,Nothing})
	@boundscheck length(Rw) == length(Kw) == length(R) ||
			error("size mismatch between position arrays")
	for i in eachindex(R)
		Rw[i], Kw[i] = pbcpos(R[i], cell)
	end
	Rw, Kw
end

struct PBCImageIterator{T<:TriclinicPBC}
	Kw::Vector3D
	cell::T

	PBCImageIterator{T}(Kw::Vector3D, cell::T) where {T<:TriclinicPBC} =
			new(Kw, cell)
end

PBCImageIterator(Kw::Vector3D, cell::T) where {T<:TriclinicPBC} =
		PBCImageIterator{T}(Kw, cell)

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
		iter.cell * (iter.Kw + pbcfirstshell[n]), n + 1
	end
end

eachimage(R::Vector3D, cell::TriclinicPBC) =
		PBCImageIterator(wrappos(inv(cell) * R, kspace), cell)

eachimage(Rw::Vector3D, Kw::Vector3D, cell::TriclinicPBC) =
		PBCImageIterator(Kw, cell)

nearestpos(R::AbstractVector{<:Real}, O::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		nearestpos(Vector3D(R), Vector3D(O), cell)

nearestpos(R::Vector3D, O::Vector3D, cell::Union{TriclinicPBC,Nothing}) =
		nearestpos(pbcpos(R, cell)..., O, cell)

function nearestpos(Rw::Vector3D, Kw::Vector3D, O::Vector3D, cell::TriclinicPBC)
	Rwmin = Rw
	d2min = sqdist(Rw, O)
	for Rwi in eachimage(Rw, Kw, cell)
		d2 = sqdist(Rwi, O)
		if d2 < d2min
			d2min = d2
			Rwmin = Rwi
		end
	end
	Rwmin
end

function nearestpos(Rw::Vector3D, Kw::Vector3D, O::Vector3D,
		cell::RhombododecahedralPBC)
	Rwmin = Rw
	d2min = sqdist(Rw, O)
	if d2min < cell.r2
		return Rwmin
	end
	for Rwi in eachimage(Rw, Kw, cell)
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

function nearestpos(Rw::Vector3D, Kw::Vector3D, O::Vector3D,
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

nearestpos(Rw::Vector3D, Kw::Vector3D, O::Vector3D, ::Nothing) = Rw

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

mindist(R1w::AbstractVector{<:Real}, K1w::AbstractVector{<:Real},
		R2w::AbstractVector{<:Real}, K2w::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		sqrt(sqmindist(R1w, K1w, R2w, K2w, cell))

sqmindist(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) =
		sqmindist(Vector3D(R1), Vector3D(R2), cell)

sqmindist(R1w::AbstractVector{<:Real}, K1w::AbstractVector{<:Real},
		R2w::AbstractVector{<:Real}, K2w::AbstractVector{<:Real},
		cell::Union{TriclinicPBC,Nothing}) = sqmindist(Vector3D(R1w),
		Vector3D(K1w), Vector3D(R2w), Vector3D(K2w), cell)

sqmindist(R1::Vector3D, R2::Vector3D, cell::Union{TriclinicPBC,Nothing}) =
		sqmindist(pbcpos(R1)..., pbcpos(R2)..., cell)

function sqmindist(R1w::Vector3D, K1w::Vector3D, R2w::Vector3D, K2w::Vector3D,
		cell::TriclinicPBC)
	d2min = sqdist(R1w, R2w)
	for R1wi in eachimage(R1w, K1w, cell)
		d2 = sqdist(R1wi, R2w)
		if d2 < d2min
			d2min = d2
		end
	end
	d2min
end

function sqmindist(R1w::Vector3D, K1w::Vector3D, R2w::Vector3D, K2w::Vector3D,
		cell::RhombododecahedralPBC)
	d2min = sqdist(R1w, R2w)
	if d2min < cell.r2
		return d2min
	end
	for R1wi in eachimage(R1w, K1w, cell)
		d2 = sqdist(R1wi, R2w)
		if d2 < cell.r2
			return d2
		elseif d2 < d2min
			d2min = d2
		end
	end
	d2min
end

function sqmindist(R1w::Vector3D, K1w::Vector3D, R2w::Vector3D, K2w::Vector3D,
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

sqmindist(R1w::Vector3D, K1w::Vector3D, R2w::Vector3D, K2w::Vector3D,
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

function unwrap!(K::AbstractVector{Vector3D}, ::Kspace,
		strategy::UnwrapByExtent)
	@boundscheck length(K) > 1 || error("cannot unwrap fewer than 2 positions")
	Kdims = dims(extent(K))
	if Kdims.x > 0.5
		for i in eachindex(K)
			x, y, z = K[i]
			x -= 0.5
			x -= floor(x)
			x += 0.5
			K[i] = Vector3D(x,y,z)
		end
	end
	if Kdims.y > 0.5
		for i in eachindex(K)
			x, y, z = K[i]
			y -= 0.5
			y -= floor(y)
			y += 0.5
			K[i] = Vector3D(x,y,z)
		end
	end
	if Kdims.z > 0.5
		for i in eachindex(K)
			x, y, z = K[i]
			z -= 0.5
			z -= floor(z)
			z += 0.5
			K[i] = Vector3D(x,y,z)
		end
	end
	K
end

struct UnwrapByGap <: UnwrapStrategy
	D::Vector3D
	buffer::Vector{Vector3D}

	function UnwrapByGap(cell::TriclinicPBC;
			d::Real = maximum(values(covalent_radii)))
		(a,b,c), angles = pbcgeometry(cell)
		new(Vector3D(d,d,d) ./ Vector3D(a,b,c), Vector3D[])
	end
end

function unwrap!(R::AbstractVector{Vector3D}, cell::TriclinicPBC,
		strategy::UnwrapByGap)
	@boundscheck length(R) > 1 || error("cannot unwrap fewer than 2 positions")
	inv(cell)(R)
	unwrap!(R, kspace, strategy)
	cell(R)
end

function unwrap!(K::AbstractVector{Vector3D}, ::Kspace,
		strategy::UnwrapByGap)
	@boundscheck length(K) > 1 || error("cannot unwrap fewer than 2 positions")
	Ks = resize!(strategy.buffer, length(K))
	Ks .= K
	sort!(Ks, by = R -> R.x)
	for i = 2:length(Ks)
		if Ks[i].x - Ks[i-1].x > strategy.D.x
			d = Ks[i].x
			for j in eachindex(K)
				x, y, z = K[j]
				x -= d
				x -= floor(x)
				x += d
				K[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	sort!(Ks, by = R -> R.y)
	for i = 2:length(Ks)
		if Ks[i].y - Ks[i-1].y > strategy.D.y
			d = Ks[i].y
			for j in eachindex(K)
				x, y, z = K[j]
				y -= d
				y -= floor(y)
				y += d
				K[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	sort!(Ks, by = R -> R.z)
	for i = 2:length(Ks)
		if Ks[i].z - Ks[i-1].z > strategy.D.z
			d = Ks[i].z
			for j in eachindex(K)
				x, y, z = K[j]
				z -= d
				z -= floor(z)
				z += d
				K[j] = Vector3D(x,y,z)
			end
			break
		end
	end
	K
end

abstract type ProximityLattice end

function proxiparams(R::AbstractVector{Vector3D}, d::Real)
	D = Vector3D(d,d,d)
	bb = extent(R)
	O = minimum(bb) - 2*D
	Dbb = dims(bb)
	nx = ceil(Int, Dbb.x / D.x) + 4
	ny = ceil(Int, Dbb.y / D.y) + 4
	nz = ceil(Int, Dbb.z / D.z) + 4
	(nx, ny, nz), O, D
end

struct NonperiodicProximityLattice <: ProximityLattice
	grid::NonperiodicGrid3D
	O::Vector3D
	D::Vector3D
end

function NonperiodicProximityLattice(R::AbstractVector{Vector3D}, d::Real)
	N, O, D = proxiparams(R, d)
	grid = NonperiodicGrid3D(N...)
	NonperiodicProximityLattice(grid, O, D)
end

function proxiparams(cell::TriclinicPBC, d::Real)
	D = proxidims(cell, d)
	nx = max(1, floor(Int, 1 / D.x))
	ny = max(1, floor(Int, 1 / D.y))
	nz = max(1, floor(Int, 1 / D.z))
	D = Vector3D(1.0/nx, 1.0/ny, 1.0/nz)
	(nx, ny, nz), D
end

function proxidims(cell::TriclinicPBC, d::Real)
	r = norm(inv(cell) * Vector3D(d, d, d))
	Vector3D(r, r, r)
end

function proxidims(cell::RhombododecahedralPBC, d::Real)
	r = √(2*d^2)
	gridcell = RhombododecahedralCell(r)
	inv(cell) * dims(gridcell)
end

proxidims(cell::OrthorhombicPBC, d::Real) = inv(cell) * Vector3D(d, d, d)

struct PeriodicProximityLattice <: ProximityLattice
	grid::PeriodicGrid3D
	D::Vector3D
end

function PeriodicProximityLattice(cell::TriclinicPBC, d::Real)
	N, D = proxiparams(cell, d)
	grid = PeriodicGrid3D(N...)
	PeriodicProximityLattice(grid, D)
end

function Base.fill!(lattice::ProximityLattice, Kw::AbstractVector{Vector3D})
	sizehint!(lattice.grid, length(Kw))
	for Kwi in Kw
		push!(lattice.grid, findcell(lattice, Kwi))
	end
	lattice
end

function Geometry.findcell(lattice::NonperiodicProximityLattice, R::Vector3D)
	x = floor(Int, (R.x - lattice.O.x) / lattice.D.x) + 1
	y = floor(Int, (R.y - lattice.O.y) / lattice.D.y) + 1
	z = floor(Int, (R.z - lattice.O.z) / lattice.D.z) + 1
	nx, ny, nz = size(lattice.grid)
	if x < 1
		x = 1
	elseif x > nx
		x = nx
	end
	if y < 1
		y = 1
	elseif y > ny
		y = ny
	end
	if z < 1
		z = 1
	elseif z > nz
		z = nz
	end
	x, y, z
end

function Geometry.findcell(lattice::PeriodicProximityLattice, Kw::Vector3D)
	x = floor(Int, Kw.x / lattice.D.x) + 1
	y = floor(Int, Kw.y / lattice.D.y) + 1
	z = floor(Int, Kw.z / lattice.D.z) + 1
	x, y, z
end

Geometry.findnear(lattice::ProximityLattice, i::Integer) =
		findnear(lattice.grid, i)

Geometry.findnear(lattice::ProximityLattice, Kw::Vector) =
		findnear(lattice.grid, findcell(lattice, Kw))

Geometry.findnear!(dest::AbstractVector{<:Integer}, lattice::ProximityLattice,
		i::Integer) = findnear!(dest, lattice.grid, i)

Geometry.findnear!(dest::AbstractVector{<:Integer}, lattice::ProximityLattice,
		Kw::Vector) = findnear!(dest, lattice.grid, findcell(lattice, Kw))

function proxilattice(R::AbstractVector{Vector3D}, ::Nothing, d::Real)
	lattice = NonperiodicProximityLattice(R, d)
	fill!(lattice, R)
end

function proxilattice(Kw::AbstractVector{Vector3D}, cell::TriclinicPBC, d::Real)
	lattice = PeriodicProximityLattice(cell, d)
	fill!(lattice, Kw)
end

end # module
