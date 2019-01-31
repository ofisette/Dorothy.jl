# Geometry-related routines, including distances and transformations

module Geometry

using LinearAlgebra
using Statistics
using ..Dorothy.Utils
using StaticArrays

export
		Vector3D, Basis3D, Matrix3D,

		sqnorm, dist, sqdist, dihedral, com, cog,

		Transformation, transform, Translation, Scaling, LinearTransformation,
		AffineTransformation, translation, scaling, rotation, superposition,
		SuperpositionStrategy, SuperposeByKabsch, FitStrategy, FitBySVD,
		fitline, fitplane, rmsd, msd,

		BoundingBox, extent, dims, center, volume, isinside,

		Grid3D, NonperiodicGrid3D, PeriodicGrid3D, findcell, findnear, findnear!

struct Vector3D <: FieldVector{3,Float64}
	x::Float64
	y::Float64
	z::Float64
end

StaticArrays.similar_type(::Type{Vector3D}, ::Type{Float64}, ::Size{(3,)}) =
		Vector3D

const Basis3D = SMatrix{3,3,Float64,9}

struct Matrix3D{T<:AbstractVector{Vector3D}} <: AbstractMatrix{Float64}
	V::T

	Matrix3D{T}(V::T) where {T<:AbstractVector{Vector3D}} = new(V)
end

Matrix3D(V::T) where {T<:AbstractVector{Vector3D}} = Matrix3D{T}(V)

Base.size(M::Matrix3D) = (3, length(M.V))

Base.getindex(M::Matrix3D, I::Vararg{Int,2}) = M.V[I[2]][I[1]]

sqnorm(V::AbstractVector{<:Real}) = sqnorm(Vector3D(V))

function sqnorm(V::Vector3D)
	d2 = V.x^2
	d2 += V.y^2
	d2 += V.z^2
end

LinearAlgebra.norm(V::Vector3D) = sqrt(sqnorm(V))

dist(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real}) =
		dist(Vector3D(R1), Vector3D(R2))

dist(R1::Vector3D, R2::Vector3D) = sqrt(sqdist(R1, R2))

sqdist(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real}) =
		sqdist(Vector3D(R1), Vector3D(R2))

function sqdist(R1::Vector3D, R2::Vector3D)
	d2 = (R2.x - R1.x)^2
	d2 += (R2.y - R1.y)^2
	d2 += (R2.z - R1.z)^2
end

Base.angle(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real}) =
		angle(Vector3D(V1), Vector3D(V2))

Base.angle(V1::Vector3D, V2::Vector3D) = acos(V1⋅V2 / (norm(V1)*norm(V2)))

Base.angle(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real},
		R3::AbstractVector{<:Real}) =
		angle(Vector3D(R1), Vector3D(R2), Vector3D(R3))

Base.angle(R1::Vector3D, R2::Vector3D, R3::Vector3D) = angle(R1-R2, R3-R2)

dihedral(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real},
		V3::AbstractVector{<:Real}) =
		dihedral(Vector3D(V1), Vector3D(V2), Vector3D(V3))

function dihedral(V1::Vector3D, V2::Vector3D, V3::Vector3D)
	V1xV2 = V1 × V2
	V2xV3 = V2 × V3
	y = (V1xV2 × V2xV3)⋅(V2/norm(V2))
	x = V1xV2 ⋅ V2xV3
	atan(y,x)
end

dihedral(R1::AbstractVector{<:Real}, R2::AbstractVector{<:Real},
		R3::AbstractVector{<:Real}, R4::AbstractVector{<:Real}) =
		dihedral(Vector3D(R1), Vector3D(R2), Vector3D(R3), Vector3D(R4))

dihedral(R1::Vector3D, R2::Vector3D, R3::Vector3D, R4::Vector3D) =
		dihedral(R2-R1, R3-R2, R4-R3)

function com(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real})
	@boundscheck begin
		length(R) > 0 || error("com of zero positions is undefined")
		length(R) == length(W) ||
				error("size mismatch between position/weight arrays")
	end
	sumx, sumy, sumz, sumW = 0.0, 0.0, 0.0, 0.0
	for i in eachindex(W)
		Wi = W[i]
		Ri = R[i]
		sumRx += Ri.x * Wi
		sumRy += Ri.y * Wi
		sumRz += Ri.z * Wi
		sumW += Wi
	end
	Vector3D(sumRx / sumW, sumRy / sumW, sumRz / sumW)
end

function cog(R::AbstractVector{Vector3D})
	@boundscheck length(R) > 0 || error("cog of zero positions is undefined")
	sumRx, sumRy, sumRz = 0.0, 0.0, 0.0
	for Ri in R
		sumRx += Ri.x
		sumRy += Ri.y
		sumRz += Ri.z
	end
	n = length(R)
	Vector3D(sumRx / n, sumRy / n, sumRz / n)
end

abstract type Transformation end

transform(T::Transformation, R::AbstractVector{<:Real}) =
		transform(T, Vector3D(R))

transform(T::Transformation, R::Vector3D) = T * R

transform(T::Transformation, R::AbstractVector{Vector3D}) = (R .= T .* R)

Base.broadcastable(T::Transformation) = Ref(T)

struct Translation <: Transformation
	V::Vector3D
end

Base.inv(T::Translation) = Translation(-T.V)

Base.:(*)(T::Translation, R::AbstractVector{<:Real}) = T * Vector3D(R)

Base.:(*)(T::Translation, R::Vector3D) = T.V + R

Base.:(*)(T2::Translation, T1::Translation) = Translation(T2.V + T1.V)

(T::Translation)(x) = transform(T, x)

struct Scaling <: Transformation
	V::Vector3D
end

Base.inv(T::Scaling) = Scaling(1.0 ./ T.V)

Base.:(*)(T::Scaling, R::AbstractVector{<:Real}) = T * Vector3D(R)

Base.:(*)(T::Scaling, R::Vector3D) = T.V .* R

Base.:(*)(T2::Scaling, T1::Scaling) = Scaling(T2.V .* T1.V)

(T::Scaling)(x) = transform(T, x)

struct LinearTransformation <:Transformation
	M::Basis3D
end

Base.inv(T::LinearTransformation) = LinearTransformation(inv(T.M))

Base.:(*)(T::LinearTransformation, R::AbstractVector{<:Real}) = T * Vector3D(R)

Base.:(*)(T::LinearTransformation, R::Vector3D) = T.M * R

Base.:(*)(T2::LinearTransformation, T1::LinearTransformation) =
		LinearTransformation(T2.M * T1.M)

Base.:(*)(T2::LinearTransformation, T1::Scaling) =
		LinearTransformation(T2.M * Diagonal(T1.V))

Base.:(*)(T2::Scaling, T1::LinearTransformation) =
		LinearTransformation(Diagonal(T2) * T1.M)

(T::LinearTransformation)(x) = transform(T, x)

struct AffineTransformation <: Transformation
	T1::LinearTransformation
	T2::Translation
end

Base.inv(T::AffineTransformation) =
		convert(AffineTransformation, inv(convert(Matrix{Float64}, T)))

Base.:(*)(T::AffineTransformation, R::AbstractVector{<:Real}) = T * Vector3D(R)

Base.:(*)(T::AffineTransformation, R::Vector3D) = T.T2 * (T.T1 * R)

Base.:(*)(T2::Translation, T1::LinearTransformation) =
		AffineTransformation(T1, T2)

Base.:(*)(T2::Translation, T1::AffineTransformation) =
		AffineTransformation(T1.T1, T2 * T1.T2)

function Base.:(*)(T2::AffineTransformation, T1::AffineTransformation)
	M2 = convert(Matrix{Float64}, T2)
	M1 = convert(Matrix{Float64}, T1)
	convert(AffineTransformation, M2 * M1)
end

(T::AffineTransformation)(x) = transform(T, x)

function Base.convert(::Type{Matrix{Float64}}, T::AffineTransformation)
	M = zeros(Float64, 4,4)
	for i = 1:3
		for j = 1:3
			M[j,i] = T.T1.M[j,i]
		end
	end
	for j = 1:3
		M[j,4] = T.T2.V[j]
	end
	M[4,4] = 1.0
	M
end

function Base.convert(::Type{AffineTransformation}, M::AbstractMatrix{<:Real})
	@boundscheck size(M) == (4,4) && M[4,:] == [0.0, 0.0, 0.0, 1.0] ||
			error("expected 4×4 affine transformation matrix")
	AffineTransformation(LinearTransformation(M[1:3,1:3]),
	Translation(M[1:3,4]))
end

Base.convert(::Type{AffineTransformation}, T::Translation) =
		AffineTransformation(LinearTransformation(Matrix{Float64}(I, 3,3)),
		Translation(copy(T.V)))

Base.convert(::Type{AffineTransformation}, T::LinearTransformation) =
		AffineTransformation(LinearTransformation(copy(T.M)),
		Translation(zeros(3)))

function Base.:(*)(T2::Transformation, T1::Transformation)
	T2p, T1p = promote(T2, T1)
	T2p * T1p
end

Base.promote_rule(::Type{<:Transformation}, ::Type{<:Transformation}) =
		AffineTransformation

translation(V::AbstractVector{<:Real}) = Translation(V)

translation(; x::Real = 0.0, y::Real = 0.0, z::Real = 0.0) =
		translation(Vector3D(x, y, z))

function scaling(V::AbstractVector{<:Real};
		O::Union{AbstractVector{<:Real},Nothing} = nothing)
	T = LinearTransformation(Diagonal(V))
	if O == nothing
		T
	else
		translation(O) * T * translation(-O)
	end
end

function scaling(c::Real; O::Union{AbstractVector{<:Real},Nothing} = nothing)
	T = scaling(Vector3D(c,c,c))
	if O == nothing
		T
	else
		translation(O) * T * translation(-O)
	end
end

function scaling(; x::Real = 1.0, y::Real = 1.0, z::Real = 1.0,
		O::Union{AbstractVector{<:Real},Nothing} = nothing)
	T = scaling(Vector3D(x,y,z))
	if O == nothing
		T
	else
		translation(O) * T * translation(-O)
	end
end

function rotation(θ::Real, V::AbstractVector{<:Real};
		O::Union{AbstractVector{<:Real},Nothing} = nothing)
	l, m, n = normalize(V)
	cosθ = cos(θ)
	sinθ = sin(θ)
	M = @SMatrix [ l*l*(1-cosθ)+cosθ  m*l*(1-cosθ)-n*sinθ  n*l*(1-cosθ)+m*sinθ ;
		          l*m*(1-cosθ)+n*sinθ  m*m*(1-cosθ)+cosθ   n*m*(1-cosθ)-l*sinθ ;
		          l*n*(1-cosθ)-m*sinθ m*n*(1-cosθ)+l*sinθ   n*n*(1-cosθ)+cosθ  ]
	T = LinearTransformation(M)
	if O == nothing
		T
	else
		translation(O) * T * translation(-O)
	end
end

function rotation(θ::Real; R1::AbstractVector{<:Real},
		R2::AbstractVector{<:Real},
		R3::Union{AbstractVector{<:Real},Nothing} = nothing)
	if R3 == nothing
		rotation(θ, R2-R1, O = R1)
	else
		V1 = R1 - R2
		V2 = R3 - R2
		rotation(θ, V1×V2, O = R2)
	end
end

function rotation(; x::Union{Real,Nothing} = nothing,
		y::Union{Real,Nothing} = nothing, z::Union{Real,Nothing} = nothing,
		O::Union{Real,Nothing} = nothing)
	if x != nothing
		if y != nothing || z != nothing
			error("expected only one of x, y or z as keyword argumen")
		end
		if O == nothing
			rotation(x, Vector3D(1.0, 0.0, 0.0))
		else
			rotation(x, Vector3D(1.0, 0.0, 0.0); O = O)
		end
	elseif y != nothing
		if x != nothing || z != nothing
			error("expected only one of x, y or z as keyword argumen")
		end
		if O == nothing
			rotation(y, Vector3D(0.0, 1.0, 0.0))
		else
			rotation(y, Vector3D(0.0, 1.0, 0.0); O = O)
		end
	elseif z != nothing
		if x != nothing || y != nothing
			error("expected only one of x, y or z as keyword argumen")
		end
		if O == nothing
			rotation(z, Vector3D(0.0, 0.0, 1.0))
		else
			rotation(z, Vector3D(0.0, 0.0, 1.0); O = O)
		end
	else
		error("expected one of x, y or z as keyword argument")
	end
end

abstract type SuperpositionStrategy end

mutable struct SuperposeByKabsch <: SuperpositionStrategy
	R::Matrix{Float64}
	Rref::Matrix{Float64}
	H::Matrix{Float64}
	T::Matrix{Float64}
	center::Bool

	SuperposeByKabsch(; center::Bool = true) =
			new(Matrix{Float64}(undef, (0,0)), Matrix{Float64}(undef, (0,0)),
			Matrix{Float64}(undef, (3,3)), Matrix{Float64}(undef, (3,3)),
			center)
end

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D}) =
		superposition(R, Rref, SuperposeByKabsch())

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D},
		W::AbstractVector{<:Real}) =
		superposition(R, Rref, W, SuperposeByKabsch())

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D},
		strategy::SuperposeByKabsch) =
		superposition(R, Rref, Repeated(1.0, length(R)), strategy)

function superposition(R::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		strategy::SuperposeByKabsch)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(Rref) == length(W) == n ||
				error("size mismatch between position/weight arrays")
	end
	if ncols(strategy.R) != n
		strategy.R = Matrix{Float64}(undef, (3,n))
		strategy.Rref = Matrix{Float64}(undef, (3,n))
	end
	for i = 1:n
		strategy.R[:,i] = R[i]
		strategy.Rref[:,i] = Rref[i]
	end
	if strategy.center
		cogR = cog(R)
		cogRref = cog(Rref)
		for i = 1:n
			strategy.R[:,i] -= cogR
			strategy.Rref[:,i] -= cogRref
		end
	end
	kabsch!(strategy.T, strategy.H, strategy.R, strategy.Rref, W)
	rotation = LinearTransformation(strategy.T)
	if strategy.center
		translation(cogRref) * rotation * translation(-cogR)
	else
		rotation
	end
end

function kabsch!(T::AbstractMatrix{<:Real}, H::AbstractMatrix{<:Real},
		R::AbstractMatrix{<:Real}, Rref::AbstractMatrix{<:Real},
		W::AbstractVector{<:Real})
	mul!(H, Rref, lmul!(Diagonal(W), R'))
	U, _, V = svd!(H)
	if sign(det(U) * det(V')) < 0
		V[:,3] .= .- V[:,3]
	end
	mul!(T, U, V')
end

abstract type FitStrategy end

mutable struct FitBySVD <:FitStrategy
	R::Matrix{Float64}
	Radj::Matrix{Float64}
	center::Bool

	FitBySVD(; center::Bool = true) = new(Matrix{Float64}(undef, (0,0)),
			Matrix{Float64}(undef, (0,0)),center)
end

fitline(R::AbstractVector{Vector3D}) = fitline(R, FitBySVD())

fitline(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real}) =
		fitline(R, W, FitBySVD())

fitline(R::AbstractVector{Vector3D}, strategy::FitBySVD) =
		fitline(R, Repeated(1.0, length(R)), strategy)

function fitline(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		strategy::FitBySVD)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(W) == n || error("size mismatch between position/weight arrays")
	end
	if ncols(strategy.R) != n
		strategy.R = Matrix{Float64}(undef, (3,n))
		strategy.Radj = Matrix{Float64}(undef, (n,3))
	end
	for i = 1:n
		strategy.R[:,i] = R[i]
	end
	if strategy.center
		cogR = cog(R)
		for i = 1:n
			strategy.R[:,i] -= cogR
		end
	end
	T = svdlinefit!(strategy.Radj, strategy.R, W)
	Vector3D(T)
end

function svdlinefit!(Radj::AbstractMatrix{<:Real}, R::AbstractMatrix{<:Real},
		W::AbstractVector{<:Real})
	mul!(Radj, Diagonal(W), R')
	T = svd!(Radj).V[:,1]
	Rn = Radj'[:,end]
	if sign(Rn ⋅ T) < 0
		T .= .- T
	end
	T
end

fitplane(R::AbstractVector{Vector3D}) = fitplane(R, FitBySVD())

fitplane(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real}) =
		fitplane(R, W, FitBySVD())

fitplane(R::AbstractVector{Vector3D}, strategy::FitBySVD) =
		fitplane(R, Repeated(1.0, length(R)), strategy)

function fitplane(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		strategy::FitBySVD)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(W) == n || error("size mismatch between position/weight arrays")
	end
	if ncols(strategy.R) != n
		strategy.R = Matrix{Float64}(undef, (3,n))
		strategy.Radj = Matrix{Float64}(undef, (n,3))
	end
	for i = 1:n
		strategy.R[:,i] = R[i]
	end
	if strategy.center
		cogR = cog(R)
		for i = 1:n
			strategy.R[:,i] -= cogR
		end
	end
	T = svdplanefit!(strategy.Radj, strategy.R, W)
	Vector3D(T)
end

function svdplanefit!(Radj::AbstractMatrix{<:Real}, R::AbstractMatrix{<:Real},
		W::AbstractVector{<:Real})
	mul!(Radj, Diagonal(W), R')
	T = svd!(Radj).V[:,3]
	R1 = Radj'[:,1]
	Rn = Radj'[:,end]
	if sign(R1 × Rn ⋅ T) < 0
		T .= .- T
	end
	T
end

rmsd(R1::AbstractVector{Vector3D}, R2::AbstractVector{Vector3D}) =
		sqrt(msd(R1, R2))

rmsd(R1::AbstractVector{Vector3D}, R2::AbstractVector{Vector3D},
		W::AbstractVector{<:Real}) = sqrt(msd(R1, R2, W))

function msd(R1::AbstractVector{Vector3D}, R2::AbstractVector{Vector3D})
	@boundscheck begin
		length(R1) > 0 || error("deviation of zero positions is undefined")
		length(R1) == length(R2) ||
				error("size mismatch between position arrays")
	end
	sd = 0.0
	for i in eachindex(R1)
		R1i = R1[i]
		R2i = R2[i]
		sd += (R2i.x - R1i.x)^2
		sd += (R2i.y - R1i.y)^2
		sd += (R2i.z - R1i.z)^2
	end
	sd / length(R1)
end

function msd(R1::AbstractVector{Vector3D}, R2::AbstractVector{Vector3D},
		W::AbstractVector{<:Real})
	@boundscheck begin
		length(R1) > 0 || error("deviation of zero positions is undefined")
		length(R1) == length(R2) == length(W) ||
				error("size mismatch between position/weight arrays")
	end
	sd = 0.0
	meanW = mean(W)
	for i in eachindex(W)
		R1i = R1[i]
		R2i = R2[i]
		Wi = W[i]
		sd += (R2i.x - R1i.x)^2 * Wi / meanW
		sd += (R2i.y - R1i.y)^2 * Wi / meanW
		sd += (R2i.z - R1i.z)^2 * Wi / meanW
	end
	sd / length(W)
end

mutable struct BoundingBox
	min::Vector3D
	max::Vector3D

	BoundingBox(min::AbstractVector{<:Real}, max::AbstractVector{<:Real}) =
			new(min, max)
end

function extent(R::AbstractVector{Vector3D})
	@boundscheck length(R) > 0 || error("extent of zero positions is undefined")
	xmin, ymin, zmin = Inf, Inf, Inf
	xmax, ymax, zmax = -Inf, -Inf, -Inf
	for Ri in R
		xmin = min(xmin, Ri.x)
		ymin = min(ymin, Ri.y)
		zmin = min(zmin, Ri.z)
		xmax = max(xmax, Ri.x)
		ymax = max(ymax, Ri.y)
		zmax = max(zmax, Ri.z)
	end
	BoundingBox(Vector3D(xmin, ymin, zmin), Vector3D(xmax, ymax, zmax))
end

extent(R::AbstractVector{<:Real}...) = extent([Vector3D(Ri) for Ri in R])

extent(bb::BoundingBox) = bb

Base.minimum(bb::BoundingBox) = bb.min

Base.maximum(bb::BoundingBox) = bb.max

dims(bb::BoundingBox) = maximum(bb) - minimum(bb)

Base.extrema(bb::BoundingBox) = minimum(bb), maximum(bb)

center(R1::Vector3D, R2::Vector3D) = R1 + (R2 - R1) / 2

center(bb::BoundingBox) = center(minimum(bb), maximum(bb))

volume(bb::BoundingBox) = prod(dims(bb))

isinside(R::AbstractVector{<:Real}, bb::BoundingBox) = isinside(Vector3D(R), bb)

isinside(R::Vector3D, bb::BoundingBox) =
		! (R.x < minimum(bb).x || R.y < minimum(bb).y || R.z < minimum(bb).z ||
		R.x > maximum(bb).x || R.y > maximum(bb).y || R.z > maximum(bb).z)

isinside(bb1::BoundingBox, bb2::BoundingBox) =
		isinside(minimum(bb1), bb2) && isinside(maximum(bb1), bb2)

abstract type Grid3D{T} end

mutable struct NonperiodicGrid3D{T} <: Grid3D{T}
	cells::Array{Vector{T},3}
	hint::Int
	N::Tuple{Int,Int,Int}
	O::Vector3D
	D::Vector3D

	function NonperiodicGrid3D{T}(bb::BoundingBox, D::Vector3D) where {T}
		g3 = new(Array{Vector{T},3}(undef, (0,0,0)), 0, (0,0,0))
		resize!(g3, bb, D)
	end
end

NonperiodicGrid3D{T}(bb::BoundingBox, D::AbstractVector{<:Real}) where {T} =
		NonperiodicGrid3D{T}(bb, Vector3D(D))

mutable struct PeriodicGrid3D{T} <: Grid3D{T}
	cells::Array{Vector{T},3}
	hint::Int
	N::Tuple{Int,Int,Int}
	O::Vector3D
	D::Vector3D

	function PeriodicGrid3D{T}(bb::BoundingBox, D::Vector3D) where {T}
		g3 = new(Array{Vector{T},3}(undef, (0,0,0)), 0, (0,0,0))
		resize!(g3, bb, D)
	end
end

PeriodicGrid3D{T}(bb::BoundingBox, D::AbstractVector{<:Real}) where {T} =
		PeriodicGrid3D{T}(bb, Vector3D(D))

Base.resize!(g3::Grid3D, bb::BoundingBox, D::AbstractVector{<:Real}) =
		resize!(g3, bb, Vector3D(D))

function Base.resize!(g3::Grid3D{T}, bb::BoundingBox, D::Vector3D) where {T}
	@boundscheck all(x -> (x > 0), D) || error("invalid cell dimension")
	(nx,ny,nz), D = gridparams(g3, bb, D)
	if any((nx,ny,nz) .> g3.N)
		g3.cells = Array{Vector{T},3}(undef, nx,ny,nz)
		for x in 1:nx
			for y in 1:ny
				for z in 1:nz
					g3.cells[x,y,z] = sizehint!(T[], g3.hint)
				end
			end
		end
	else
		empty!(g3)
	end
	g3.N = (nx,ny,nz)
	g3.O = minimum(bb)
	g3.D = D
	g3
end

function gridparams(g3::NonperiodicGrid3D, bb::BoundingBox, D::Vector3D)
	N = floor.(Int, dims(bb) ./ D) .+ 1
	N, D
end

function gridparams(g3::PeriodicGrid3D, bb::BoundingBox, D::Vector3D)
	N = max.(1, floor.(Int, dims(bb) ./ D))
	D = dims(bb) ./ N
	N, D
end

function Base.sizehint!(g3::Grid3D, n::Integer)
	g3.hint = n
	for i in eachindex(g3.cells)
		sizehint!(g3.cells[i], n)
	end
	g3
end

function Base.empty!(g3::Grid3D)
	(nx,ny,nz) = g3.N
	for x in 1:nx
		for y in 1:ny
			for z in 1:nz
				empty!(g3.cells[x,y,z])
			end
		end
	end
	g3
end

Base.push!(g3::Grid3D{T}, R::AbstractVector{<:Real}, v::T) where {T} =
		push!(g3, Vector3D(R), v)

Base.push!(g3::Grid3D{T}, R::Vector3D, v::T) where {T} =
		push!(g3, findcell(g3, R), v)

function Base.push!(g3::Grid3D{T}, (x,y,z)::Tuple{Integer,Integer,Integer},
		v::T) where {T}
	#=
	@boundscheck begin
		(nx,ny,nz) = g3.N
		(1 <= x <= nx) && (1 <= y <= ny) && (1 <= z <= nz) ||
				error("cannot push value to a cell outside the grid")
	end
	=#
	push!(g3.cells[x,y,z], v)
	g3
end

findcell(g3::Grid3D, R::AbstractVector{<:Real}) = findcell(g3, Vector3D(R))

function findcell(g3::NonperiodicGrid3D, R::Vector3D)
	x = floor(Int, (R.x - g3.O.x) / g3.D.x) + 1
	y = floor(Int, (R.y - g3.O.y) / g3.D.y) + 1
	z = floor(Int, (R.z - g3.O.z) / g3.D.z) + 1
	x, y, z
end

function findcell(g3::PeriodicGrid3D, R::Vector3D)
	x = floor(Int, (R.x - g3.O.x) / g3.D.x) + 1
	y = floor(Int, (R.y - g3.O.y) / g3.D.y) + 1
	z = floor(Int, (R.z - g3.O.z) / g3.D.z) + 1
	(nx,ny,nz) = g3.N
	if x > nx
		x = 1
	end
	if y > ny
		y = 1
	end
	if z > nz
		z = 1
	end
	x, y, z
end

findnear(g3::Grid3D{T}, R::AbstractVector{<:Real}) where {T} =
		findnear(T[], g3, R)

findnear(g3::Grid3D{T}, (x,y,z)::Tuple{Integer,Integer,Integer}) where {T} =
		findnear!(T[], g3, (x,y,z))

findnear!(dest::AbstractVector{T}, g3::Grid3D{T},
		R::AbstractVector{<:Real}) where {T} = findnear!(dest, g3, Vector3D(R))

findnear!(dest::AbstractVector{T}, g3::Grid3D{T}, R::Vector3D) where {T} =
		findnear!(dest, g3, findcell(R))

function findnear!(dest::AbstractVector{T}, g3::NonperiodicGrid3D{T},
		(x,y,z)::Tuple{Integer,Integer,Integer}) where {T}
	empty!(dest)
	(nx,ny,nz) = g3.N
	for xi = x-1:x+1
		for yi = y-1:y+1
			for zi = z-1:z+1
				if 1 <= xi <= nx && 1 <= yi <= ny && 1 <= zi <= nz
					for v in g3.cells[xi,yi,zi]
						push!(dest, v)
					end
				end
			end
		end
	end
	dest
end

function findnear!(dest::AbstractVector{T}, g3::PeriodicGrid3D{T},
		(x,y,z)::Tuple{Integer,Integer,Integer}) where {T}
	empty!(dest)
	(nx,ny,nz) = g3.N
	for xi = x-1:x+1
		for yi = y-1:y+1
			for zi = z-1:z+1
				if xi == 0
					xi = nx
				elseif xi == nx + 1
					xi = 1
				end
				if yi == 0
					yi = ny
				elseif yi == ny + 1
					yi = 1
				end
				if zi == 0
					zi = nz
				elseif zi == nz + 1
					zi = 1
				end
				for v in g3.cells[xi,yi,zi]
					push!(dest, v)
				end
			end
		end
	end
	dest
end

end # module
