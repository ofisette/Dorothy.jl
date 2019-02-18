# Geometry-related routines, including distances and transformations

module Geometry

using LinearAlgebra
using Statistics
using ..Dorothy.Utils
using StaticArrays

export
		Vector3D, Basis3D,

		sqnorm, dist, sqdist, dihedral, com, cog,

		Transformation, transform, Translation, Scaling, LinearTransformation,
		AffineTransformation, translation, scaling, rotation, superposition,
		SuperpositionBuffer, FitBuffer, fitline, fitplane, rmsd, msd,

		BoundingBox, extent, dims, center, volume, isinside

struct Vector3D <: FieldVector{3,Float64}
	x::Float64
	y::Float64
	z::Float64
end

StaticArrays.similar_type(::Type{Vector3D}, ::Type{Float64}, ::Size{(3,)}) =
		Vector3D

const Basis3D = SMatrix{3,3,Float64,9}

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
	M = @SMatrix [l*l*(1-cosθ)+cosθ   m*l*(1-cosθ)-n*sinθ n*l*(1-cosθ)+m*sinθ;
		          l*m*(1-cosθ)+n*sinθ m*m*(1-cosθ)+cosθ   n*m*(1-cosθ)-l*sinθ;
		          l*n*(1-cosθ)-m*sinθ m*n*(1-cosθ)+l*sinθ n*n*(1-cosθ)+cosθ  ]
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
			error("expected only one of x, y or z as keyword argument")
		end
		if O == nothing
			rotation(x, Vector3D(1.0, 0.0, 0.0))
		else
			rotation(x, Vector3D(1.0, 0.0, 0.0); O = O)
		end
	elseif y != nothing
		if x != nothing || z != nothing
			error("expected only one of x, y or z as keyword argument")
		end
		if O == nothing
			rotation(y, Vector3D(0.0, 1.0, 0.0))
		else
			rotation(y, Vector3D(0.0, 1.0, 0.0); O = O)
		end
	elseif z != nothing
		if x != nothing || y != nothing
			error("expected only one of x, y or z as keyword argument")
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

mutable struct SuperpositionBuffer
	R::Matrix{Float64}
	Rref::Matrix{Float64}
	H::Matrix{Float64}
	T::Matrix{Float64}
	center::Bool

	SuperpositionBuffer(; center::Bool = true) =
			new(Matrix{Float64}(undef, (0,0)), Matrix{Float64}(undef, (0,0)),
			Matrix{Float64}(undef, (3,3)), Matrix{Float64}(undef, (3,3)),
			center)
end

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D}) =
		superposition(R, Rref, SuperpositionBuffer())

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D},
		W::AbstractVector{<:Real}) =
		superposition(R, Rref, W, SuperpositionBuffer())

superposition(R::AbstractVector{Vector3D}, Rref::AbstractVector{Vector3D},
		buffer::SuperpositionBuffer) =
		superposition(R, Rref, Repeated(1.0, length(R)), buffer)

function superposition(R::AbstractVector{Vector3D},
		Rref::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		buffer::SuperpositionBuffer)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(Rref) == length(W) == n ||
				error("size mismatch between position/weight arrays")
	end
	if ncols(buffer.R) != n
		buffer.R = Matrix{Float64}(undef, (3,n))
		buffer.Rref = Matrix{Float64}(undef, (3,n))
	end
	for i = 1:n
		buffer.R[:,i] = R[i]
		buffer.Rref[:,i] = Rref[i]
	end
	if buffer.center
		cogR = cog(R)
		cogRref = cog(Rref)
		for i = 1:n
			buffer.R[:,i] -= cogR
			buffer.Rref[:,i] -= cogRref
		end
	end
	kabsch!(buffer.T, buffer.H, buffer.R, buffer.Rref, W)
	rotation = LinearTransformation(buffer.T)
	if buffer.center
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

mutable struct FitBuffer
	R::Matrix{Float64}
	Radj::Matrix{Float64}
	center::Bool

	FitBuffer(; center::Bool = true) = new(Matrix{Float64}(undef, (0,0)),
			Matrix{Float64}(undef, (0,0)),center)
end

fitline(R::AbstractVector{Vector3D}) = fitline(R, FitBuffer())

fitline(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real}) =
		fitline(R, W, FitBuffer())

fitline(R::AbstractVector{Vector3D}, buffer::FitBuffer) =
		fitline(R, Repeated(1.0, length(R)), buffer)

function fitline(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		buffer::FitBuffer)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(W) == n || error("size mismatch between position/weight arrays")
	end
	if ncols(buffer.R) != n
		buffer.R = Matrix{Float64}(undef, (3,n))
		buffer.Radj = Matrix{Float64}(undef, (n,3))
	end
	for i = 1:n
		buffer.R[:,i] = R[i]
	end
	if buffer.center
		cogR = cog(R)
		for i = 1:n
			buffer.R[:,i] -= cogR
		end
	end
	T = svdlinefit!(buffer.Radj, buffer.R, W)
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

fitplane(R::AbstractVector{Vector3D}) = fitplane(R, FitBuffer())

fitplane(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real}) =
		fitplane(R, W, FitBuffer())

fitplane(R::AbstractVector{Vector3D}, buffer::FitBuffer) =
		fitplane(R, Repeated(1.0, length(R)), buffer)

function fitplane(R::AbstractVector{Vector3D}, W::AbstractVector{<:Real},
		buffer::FitBuffer)
	n = length(R)
	@boundscheck begin
		n > 0 || error("fit of zero positions is undefined")
		length(W) == n || error("size mismatch between position/weight arrays")
	end
	if ncols(buffer.R) != n
		buffer.R = Matrix{Float64}(undef, (3,n))
		buffer.Radj = Matrix{Float64}(undef, (n,3))
	end
	for i = 1:n
		buffer.R[:,i] = R[i]
	end
	if buffer.center
		cogR = cog(R)
		for i = 1:n
			buffer.R[:,i] -= cogR
		end
	end
	T = svdplanefit!(buffer.Radj, buffer.R, W)
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

struct BoundingBox
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

end # module
