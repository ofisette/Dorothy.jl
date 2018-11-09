# Geometry-related routines, including distances and transformations

module Geometry

using LinearAlgebra
using Statistics
using ..Dorothy.Utils

export
        sqnorm, dist, sqdist, dihedral, com, com!, cog, cog!, msd, rmsd,

        Transformation, translation, rotation, scaling,

        superposition, kabsch, kabsch!, fitline, fitline!, fitplane, fitplane!,

        BoundingBox, extent, extent!, dim, center,

        PBCCell, pbccell, pbcbox, volume, iscubic, ishexagonal,
		istetragonal, isorthorhombic, ismonoclinic, isinside, isinside!,

		wrappos!, wrapcog!, unwrappedmaxdim, unwrapbydim!, mindist, sqmindist,
        compactpos!, compactcog!, sqmindist!, mintrans!,

        OrthoGrid3, KspaceGrid3, findcell, findnear, findnear!

function sqnorm(V::AbstractVector{<:Real})
    d2 = 0.0
    for i in eachindex(V)
        d2 += V[i]^2
    end
    d2
end

dist(A::AbstractVector{<:Real}, B::AbstractVector{<:Real}) =
        sqrt(sqdist(A, B))

function sqdist(A::AbstractVector{<:Real}, B::AbstractVector{<:Real})
    d2 = 0.0
    for i = eachindex(A)
        d2 += (A[i] - B[i])^2
    end
    d2
end

Base.angle(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real}) =
        acos(V1⋅V2 / (norm(V1)*norm(V2)))

Base.angle(A::AbstractVector{<:Real}, B::AbstractVector{<:Real},
        C::AbstractVector{<:Real}) = angle(A-B, C-B)

function dihedral(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real},
        V3::AbstractVector{<:Real})
    V1xV2 = V1 × V2
    V2xV3 = V2 × V3
    y = (V1xV2 × V2xV3)⋅(V2/norm(V2))
    x = V1xV2 ⋅ V2xV3
    atan(y,x)
end

dihedral(A::AbstractVector{<:Real}, B::AbstractVector{<:Real},
        C::AbstractVector{<:Real}, D::AbstractVector{<:Real}) =
        dihedral(B-A, C-B, D-C)

com(R::AbstractMatrix{<:Real}, W::AbstractVector{<:Real}) =
        com!(similar(R, nrows(R)), R, W)

function com!(dest::AbstractVector{<:Real}, R::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real})
    dest .= 0.0
    sumW = 0.0
    nrows, ncols = size(R)
    for i = 1:ncols
        Wi = W[i]
        sumW += Wi
        for j = 1:nrows
            dest[j] += R[j,i] * Wi
        end
    end
    dest ./= sumW
end

cog(R::AbstractMatrix{<:Real}) = cog!(similar(R, nrows(R)), R)

cog!(dest::AbstractVector{<:Real}, R::AbstractMatrix{<:Real}) =
        com!(dest, R, ScalarArray(1.0, ncols(R)))

rmsd(R1::AbstractMatrix{<:Real}, R2::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real} = ScalarArray(1.0, ncols(R1))) =
        sqrt(msd(R1, R2, W))

function msd(R1::AbstractMatrix{<:Real}, R2::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real} = ScalarArray(1.0, ncols(R1)))
    nrows, ncols = size(R1)
    @boundscheck begin
        size(R2) == (nrows, ncols) ||
                error("dimension mismatch between position matrices")
        length(W) == ncols ||
                error("dimension mismatch between position and weight arrays")
    end
    sd = 0.0
    meanW = mean(W)
    for i = 1:ncols
        for j = 1:nrows
            d = R2[j,i] - R1[j,i]
            sd += d^2 * W[i] / meanW
        end
    end
    sd / ncols
end

abstract type Transformation end

struct Translation <: Transformation
    V::Vector{Float64}

    function Translation(V::AbstractVector{<:Real})
        @boundscheck length(V) == 3 || error("expected 3-element vector")
        new(V)
    end
end

(T::Translation)(R::AbstractArray{<:Real}) = T.V .+ R

(T2::Translation)(T1::Transformation) = T2 * T1

(T2::Translation)(T1::Translation) = Translation(T2.V + T1.V)

LinearAlgebra.mul!(R′::AbstractArray{<:Real}, T::Translation,
        R::AbstractArray{<:Real}) = (R′ .= T.V .+ R)

LinearAlgebra.lmul!(T::Translation, R::AbstractArray{<:Real}) = (R .+= T.V)

struct LinearTransformation <:Transformation
    M::Matrix{Float64}

    function LinearTransformation(M::AbstractMatrix{<:Real})
        @boundscheck size(M) == (3,3) || error("expected 3×3 matrix")
        new(M)
    end
end

(T::LinearTransformation)(R::AbstractArray{<:Real}) = T.M * R

(T2::LinearTransformation)(T1::Transformation) = T2 * T1

(T2::LinearTransformation)(T1::LinearTransformation) =
        LinearTransformation(T2.M * T1.M)

LinearAlgebra.mul!(R′::AbstractArray{<:Real}, T::LinearTransformation,
        R::AbstractArray{<:Real}) = mul!(R′, T.M, R)

struct AffineTransformation <: Transformation
    T1::LinearTransformation
    T2::Translation

    AffineTransformation(T1::LinearTransformation, T2::Translation) =
            new(T1, T2)
end

function Base.convert(::Type{Matrix{T2}}, T::AffineTransformation) where
        {T2<:Real}
    M = zeros(T2, 4,4)
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

Base.convert(::Type{AffineTransformation}, M::AbstractMatrix{<:Real}) =
        AffineTransformation(LinearTransformation(M[1:3,1:3]),
        Translation(M[1:3,4]))

Base.convert(::Type{AffineTransformation}, T::Translation) =
        AffineTransformation(LinearTransformation(Matrix{Float64}(I, 3,3)),
        Translation(copy(T.V)))

Base.convert(::Type{AffineTransformation}, T::LinearTransformation) =
        AffineTransformation(LinearTransformation(copy(T.M)),
        Translation(zeros(3)))

(T::AffineTransformation)(R::AbstractArray{<:Real}) = T.T2(T.T1(R))

(T2::AffineTransformation)(T1::Transformation) = T2 * T1

function (T2::AffineTransformation)(T1::AffineTransformation)
    M2 = convert(Matrix{Float64}, T2)
    M1 = convert(Matrix{Float64}, T1)
    convert(AffineTransformation, M2*M1)
end

function LinearAlgebra.mul!(R′::AbstractArray{<:Real}, T::AffineTransformation,
        R::AbstractArray{<:Real})
    mul!(R′, T.T1, R)
    mul!(R′, T.T2, R)
end

Base.:(*)(T::Transformation, A::AbstractArray{<:Real}) = T(A)

function Base.:(*)(T2::Transformation, T1::Transformation)
    T2′, T1′ = promote(T2, T1)
    T2′(T1′)
end

Base.promote_rule(::Type{<:Transformation}, ::Type{<:Transformation}) =
        AffineTransformation

translation(T::AbstractVector{<:Real}) = Translation(copy(T))

translation(x::Real, y::Real, z::Real) = Translation([x,y,z])

translation(; x::Real = 0.0, y::Real = 0.0, z::Real = 0.0) =
        Translation([x,y,z])

function rotation(θ::Real, V::AbstractVector{<:Real})
    l, m, n = normalize(V)
    cosθ = cos(θ)
    sinθ = sin(θ)
    M = [  l*l*(1-cosθ)+cosθ   m*l*(1-cosθ)-n*sinθ  n*l*(1-cosθ)+m*sinθ ;
          l*m*(1-cosθ)+n*sinθ   m*m*(1-cosθ)+cosθ   n*m*(1-cosθ)-l*sinθ ;
          l*n*(1-cosθ)-m*sinθ  m*n*(1-cosθ)+l*sinθ   n*n*(1-cosθ)+cosθ  ]
    LinearTransformation(M)
end

rotation(θ::Real, V::AbstractVector{<:Real}, O::AbstractVector{<:Real}) =
        Translation(O) * rotation(θ, V) * Translation(-O)

@inline rotation(θ::Real, axis::Symbol) = rotation(θ, Val(axis))

@inline rotation(θ::Real, ::Val{:x}) = rotation(θ, [1.0,0.0,0.0])

@inline rotation(θ::Real, ::Val{:y}) = rotation(θ, [0.0,1.0,0.0])

@inline rotation(θ::Real, ::Val{:z}) = rotation(θ, [0.0,0.0,1.0])

scaling(V::AbstractVector) = LinearTransformation(collect(Diagonal(V)))

scaling(x::Real, y::Real, z::Real) = scaling([x,y,z])

scaling(c::Real) = scaling(c, c, c)

@inline function superposition(R::AbstractMatrix{<:Real},
        Rref::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real} = ScalarArray(1.0, ncols(R));
        alg::Symbol = :kabsch, kwargs...)
    @boundscheck begin
        size(R) == size(Rref) || error("incompatible position matrices")
        length(W) == ncols(R) || error("incompatible weight vector")
        nrows(R) == nrows(Rref) == 3 || error("expected 3-coordinate positions")
    end
    superposition(Val(alg), R, Rref, W; kwargs...)
end

@inline superposition(::Val{alg}, R::AbstractMatrix{<:Real},
        Rref::AbstractMatrix{<:Real}, W::AbstractVector{<:Real};
        kwargs...) where {alg} =
        error("unknown superposition algorithm: $(alg)")

@inline superposition(::Val{:kabsch}, R::AbstractMatrix{<:Real},
        Rref::AbstractMatrix{<:Real}, W::AbstractVector{<:Real}; kwargs...) =
        kabsch(R, Rref, W; kwargs...)

function kabsch(R::AbstractMatrix{<:Real}, Rref::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real})
    H = similar(R, 3,3)
    comR = similar(R, 3)
    comRref = similar(Rref, 3)
    Rc = similar(R)
    Rrefc = similar(Rref)
    kabsch!(H, comR, comRref, Rc, Rrefc, R, Rref, W)
end

function kabsch!(H::AbstractMatrix{<:Real}, comR::AbstractVector{<:Real},
        comRref::AbstractVector{<:Real}, Rc::AbstractMatrix{<:Real},
        Rrefc::AbstractMatrix{<:Real}, R::AbstractMatrix{<:Real},
        Rref::AbstractMatrix{<:Real}, W::AbstractVector{<:Real})
    Rc .= R
    Rrefc .= Rref
    com!(comR, Rc, W)
    com!(comRref, Rrefc, W)
    Rc .-= comR
    Rrefc .-= comRref
    mul!(H, Rrefc, lmul!(Diagonal(W), Rc'))
    U, _, V = svd(H)
    if sign(det(V') * det(U)) < 0
        V'[3,:] .= -V'[3,:]
    end
    Translation(comRref) * LinearTransformation(U * V') * Translation(-comR)
end

function fitline(R::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real} = ScalarArray(1.0, ncols(R)))
    @boundscheck begin
        nrows, ncols = size(R)
        nrows == 3 || error("expected 3-coordinate positions")
        ncols >= 2 || error("line fitting requires at least 2 points")
        length(W) == ncols || error("incompatible weight vector")
    end
    V = similar(R, 3)
    comR = similar(R, 3)
    Rc = similar(R)
    fitline!(V, comR, Rc, R, W)
end

function fitline!(V::AbstractVector{<:Real}, comR::AbstractVector{<:Real},
        Rc::AbstractMatrix{<:Real}, R::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real})
    Rc .= R
    com!(comR, Rc, W)
    Rc .-= comR
    lmul!(Diagonal(W), Rc')
    V .= svd(Rc').V[:,1]
    Rcn = @view Rc[:,end]
    if sign(Rcn ⋅ V) < 0
        V .= .- V
    end
    V, comR
end

function fitplane(R::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real} = ScalarArray(1.0, ncols(R)))
    @boundscheck begin
        nrows, ncols = size(R)
        nrows == 3 || error("expected 3-coordinate positions")
        ncols >= 3 || error("plane fitting requires at least 3 points")
        length(W) == ncols || error("incompatible weight vector")
    end
    V = similar(R, 3)
    comR = similar(R, 3)
    Rc = similar(R)
    fitplane!(V, comR, Rc, R, W)
end

function fitplane!(V::AbstractVector{<:Real}, comR::AbstractVector{<:Real},
        Rc::AbstractMatrix{<:Real}, R::AbstractMatrix{<:Real},
        W::AbstractVector{<:Real})
    Rc .= R
    com!(comR, Rc, W)
    Rc .-= comR
    lmul!(Diagonal(W), Rc')
    V .= svd(Rc').V[:,3]
    Rc1 = @view Rc[:,1]
    Rcn = @view Rc[:,end]
    if sign(Rc1 × Rcn ⋅ V) < 0
        V .= .- V
    end
    V, comR
end

mutable struct BoundingBox{T<:Real}
    minimum::Vector{T}
    maximum::Vector{T}

    BoundingBox{T}(minimum::AbstractVector{T},
            maximum::AbstractVector{T}) where {T<:Real} = new(minimum, maximum)
end

BoundingBox(minimum::AbstractVector{T}, maximum::AbstractVector{T}) where
        {T<:Real} = BoundingBox{T}(minimum, maximum)

BoundingBox{T}(::UndefInitializer, n::Integer) where {T<:Real} =
        BoundingBox(Vector{T}(undef, n), Vector{T}(undef, n))

extent(R::AbstractMatrix{T}) where {T<:Real} =
        extent!(BoundingBox{T}(undef, nrows(R)), R)

function extent!(bb::BoundingBox{<:Real}, R::AbstractMatrix{<:Real})
    minimum = bb.minimum .= R[:,1]
    maximum = bb.maximum .= R[:,1]
    nrows, ncols = size(R)
    for j = 2:ncols
        for i = 1:nrows
            minimum[i] = min(minimum[i], R[i,j])
            maximum[i] = max(maximum[i], R[i,j])
        end
    end
    bb
end

dim(bb::BoundingBox) = bb.maximum - bb.minimum

Base.minimum(bb::BoundingBox) = copy(bb.minimum)

Base.maximum(bb::BoundingBox) = copy(bb.maximum)

function center(bb::BoundingBox)
    A = dim(bb)
    A ./= 2
    A .+= bb.minimum
end

function isinside(A::AbstractVector{<:Real}, bb::BoundingBox)
	for i in eachindex(A)
		if A < bb.minimum[i] || A > bb.maximum[i]
			return false
		end
	end
	true
end

function isinside(bb1::BoundingBox, bb2::BoundingBox)
	for i in eachindex(bb2.minimum)
		if bb1.minimum[i] < bb2.minimum[i] || bb1.maximum[i] > bb2.maximum[i]
			return false
		end
	end
	true
end

struct PBCCell{T<:Real} <: AbstractMatrix{T}
	M::Matrix{T}
	inv::Matrix{T}

	@inline function PBCCell{T}(M::AbstractMatrix{T}) where {T<:Real}
		@boundscheck size(M) == (3,3) || error("expected 3×3 matrix")
		new(M, inv(M))
	end
end

PBCCell(M::AbstractMatrix{T}) where {T<:Real} = PBCCell{T}(M)

Base.convert(::Type{PBCCell}, cell::PBCCell) = cell

Base.convert(::Type{PBCCell}, M::AbstractMatrix{<:Real}) = PBCCell(M)

Base.size(cell::PBCCell) = (3,3)

Base.getindex(cell::PBCCell, i::Int) = cell.M[i]

Base.IndexStyle(::Type{<:PBCCell}) = IndexLinear()

Base.inv(cell::PBCCell) = copy(cell.inv)

pbccell(xyz::Real) = pbccell([xyz, xyz, xyz], [τ/4, τ/4, τ/4])

pbccell(xy::Real, z::Real) = pbccell([xy, xy, z], [τ/4, τ/4, τ/4])

pbccell(x::Real, y::Real, z::Real) = pbccell([x, y, z], [τ/4, τ/4, τ/4])

function pbccell((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real})
	Ax, Ay, Az = a, 0.0, 0.0 # A lies along the positive X axis
	Bx, By, Bz = b*cos(α), b*sin(α), 0.0 # B lies in the X-Y plane
	Cx = c*cos(β)
	Cy = (b*c*cos(γ) - Bx*Cx) / By
	Cz = √(c^2 - Cx^2 - Cy^2)
	PBCCell([[Ax,Ay,Az] [Bx,By,Bz] [Cx,Cy,Cz]])
end

pbccell(M::AbstractMatrix{<:Real}) = PBCCell(collect(Float64, M))

pbccell(cell::PBCCell) = PBCCell(cell.M)

pbccell(::UniformScaling) = PBCCell(Matrix{Float64}(I, 3,3))

function pbcbox(cell::PBCCell)
    A, B, C = [cell[1:3,i] for i = 1:3]
    a, b, c = [norm(V) for V in (A, B, C)]
    α, β, γ = acos(B⋅C / (b*c)), acos(A⋅C / (a*c)), acos(A⋅B / (a*b))
    (lengths=[a,b,c], angles=[α,β,γ])
end

function volume(cell::PBCCell)
	A, B, C = [cell[1:3,i] for i = 1:3]
	A⋅(B×C)
end

function volume((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real})
	cosα, cosβ, cosγ = cos(α), cos(β), cos(γ)
	a*b*c * √(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
end

function iscubic((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real})
	a == b == c || return false
	α == β == γ == τ / 4 || return false
	true
end

function ishexagonal((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real})
	a == b || return false
	α == β == τ / 4 || return false
	γ == τ / 3 || return false
	true
end

function istetragonal((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real})
	a == b || return false
	α == β == γ == τ / 4 || return false
	true
end

isorthorhombic((a,b,c)::AbstractVector{<:Real},
		(α,β,γ)::AbstractVector{<:Real}) = (α == β == γ == τ / 4)

ismonoclinic((a,b,c)::AbstractVector{<:Real}, (α,β,γ)::AbstractVector{<:Real}) =
		(α == τ / 4)

function isinside(A::AbstractVector{<:Real}, cell::PBCCell)
    Ak = similar(A)
    isinside!(Ak, A, cell)
end

function isinside!(Ak::AbstractVector{<:Real}, A::AbstractVector{<:Real},
        cell::PBCCell)
    mul!(Ak, cell.inv, A)
    isinside(Ak)
end

function isinside(Ak::AbstractVector{<:Real})
    for i in eachindex(Ak)
        if ! (0.0 <= Ak[i] < 1.0)
            return false
        end
    end
    true
end

function wrappos!(R::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
		cell::PBCCell)
    Ak = similar(R, 3)
    wrappos!(R, Ak, cell)
end

function wrappos!(R::AbstractMatrix{<:Real}, Rik::AbstractVector{<:Real},
        cell::PBCCell)
    for i = 1:ncols(R)
        Ri = @view R[:,i]
		wrappos!(Ri, Rik, cell)
    end
    R
end

function wrappos!(A::AbstractVector{<:Real}, Ak::AbstractVector{<:Real},
        cell::PBCCell)
    mul!(Ak, cell.inv, A)
    wrappos!(Ak)
    mul!(A, cell, Ak)
end

wrappos!(Rk::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}) =
        (Rk .-= floor.(Rk))

function wrapcog!(R::AbstractMatrix{T}, cell::PBCCell) where {T<:Real}
    cogR = similar(R, 3)
    cogRk = similar(R, 3)
    wrapcog!(R, cogR, cogRk, cell)
end

function wrapcog!(R::AbstractMatrix{<:Real}, cogR::AbstractVector{<:Real},
        cogRk::AbstractVector{<:Real}, cell::PBCCell)
    cog!(cogR, R)
    R .-= cogR
    wrappos!(cogR, cogRk, cell)
    R .+= cogR
end

unwrappedmaxdim(cell::PBCCell, dmax::Real = 2.5) =
    	cell.inv * [dmax, dmax, dmax]

function unwrapbydim!(R::AbstractMatrix{<:Real}, cell::PBCCell,
        Dkmax::AbstractVector{<:Real})
	Rk = similar(R)
    Rksi = similar(R, ncols(R))
	unwrapbydim!(R, Rk, Rksi, cell, Dkmax)
end

function unwrapbydim!(R::AbstractMatrix{<:Real}, Rk::AbstractMatrix{<:Real},
        Rksi::AbstractVector{<:Real}, cell::PBCCell,
        Dkmax::AbstractVector{<:Real})
	mul!(Rk, cell.inv, R)
	unwrapbydim!(Rk, Rksi, Dkmax)
	mul!(R, cell, Rk)
end

function unwrapbydim!(Rk::AbstractMatrix{<:Real}, Rksi::AbstractVector{<:Real},
		Dkmax::AbstractVector{<:Real})
    for i in eachindex(Dkmax)
        Rki = @view Rk[i,:]
        minRki, maxRki = extrema(Rki)
        if maxRki - minRki > 1.0 - Dkmax[i]
            Rki .-= 0.5
            wrappos!(Rki)
            Rki .+= 0.5
            minRki, maxRki = extrema(Rki)
            if maxRki - minRki > 1.0 - Dkmax[i]
                wrappos!(Rki)
                Rksi .= Rki
                sort!(Rksi)
                k = 0
                for j = 1:length(Rksi)-1
                    if Rksi[j+1] - Rksi[j] > Dkmax[i]
                        k = j
                        break
                    end
                end
                if k == 0
                    error("failed to unwrap from dimensions")
                else
                    t = Rksi[k] + Dkmax[i] / 2
                    Rki .-= t
                    wrappos!(Rki)
                    Rki .+= t
                end
            end
        end
    end
end

mindist(A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, cell::PBCCell) =
		sqrt(sqmindist(A, B, cell))

function sqmindist(A::AbstractVector{<:Real}, B::AbstractVector{<:Real},
		cell::PBCCell)
	Ak = similar(A, 3)
	Bk = similar(A, 3)
	Tk = similar(A, 3)
	T = similar(A, 3)
	sqmindist!(Tk, Ak, Bk, T, A, B, cell)
end

function sqmindist!(Tk::AbstractVector{<:Real}, Ak::AbstractVector{<:Real},
		Bk::AbstractVector{<:Real}, T::AbstractVector{<:Real},
		A::AbstractVector{<:Real}, B::AbstractVector{<:Real},
		cell::PBCCell)
	mul!(Ak, cell.inv, A)
	mul!(Bk, cell.inv, B)
	mintrans!(Tk, Ak, Bk)
	mul!(T, cell, Tk)
	sqnorm(T)
end

function mintrans!(Tk::AbstractVector{<:Real}, Ak::AbstractVector{<:Real},
		Bk::AbstractVector{<:Real})
	Tk .= abs.(Bk .- Ak)
	Tk .= min.(Tk, 1.0 .- Tk)
end

function compactpos!(R::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
		cell::PBCCell)
    Ak = similar(R, 3)
    compactpos!(R, Ak, cell)
end

function compactpos!(R::AbstractMatrix{<:Real}, Rik::AbstractVector{<:Real},
        cell::PBCCell)
    for i = 1:ncols(R)
        Ri = @view R[:,i]
        compactpos!(Ri, Rik, cell)
    end
    R
end

function compactpos!(A::AbstractVector{<:Real}, Ak::AbstractVector{<:Real},
        cell::PBCCell)
    mul!(Ak, cell.inv, A)
    compactpos!(Ak)
    mul!(A, cell, Ak)
end

function compactpos!(Rk::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}})
    # wrappos!(Rk)
	Rk .+= 0.5
	wrappos!(Rk)
	Rk .-= 0.5
end

function compactcog!(R::AbstractMatrix{<:Real}, cell::PBCCell)
    Rc = similar(R, 3)
    Rck = similar(R, 3)
    compactcog!(R, Rc, Rck, cell)
end

function compactcog!(R::AbstractMatrix{<:Real}, Rc::AbstractVector{<:Real},
        Rck::AbstractVector{<:Real}, cell::PBCCell)
    cog!(Rc, R)
    R .-= Rc
    compactcog!(Rc, Rck, cell)
    R .+= Rc
end

abstract type Grid3{T} <: AbstractArray{T,3}
    #=
    struct
        cells::Array{Vector{T},3}
    end
    =#
end

struct OrthoGrid3{T} <: Grid3{T}
	cells::Array{Vector{T},3}
    Ox::Float64
    Oy::Float64
    Oz::Float64
    d::Float64

	function OrthoGrid3{T}(bb::BoundingBox, d::Real) where {T}
		@boundscheck d > 0.0 || error("invalid cell dimension")
        O = minimum(bb)
        S = dim(bb)
        N = floor.(Int, S ./ d) .+ 1
        cells = Array{Vector{T},3}(undef, N...)
        for i in eachindex(cells)
            cells[i] = T[]
        end
        new(cells, O[1], O[2], O[3], d)
	end
end

struct KspaceGrid3{T} <: Grid3{T}
	cells::Array{Vector{T},3}
    Dx::Float64
    Dy::Float64
    Dz::Float64

	function KspaceGrid3{T}(D::AbstractVector{<:Real}) where {T}
		@boundscheck begin
			for d in D
                d > 0.0 || error("invalid cell dimension")
            end
		end
        N = floor.(Int, 1.0 ./ D)
        D = 1.0 ./ N
        cells = Array{Vector{T},3}(undef, N...)
        for i in eachindex(cells)
            cells[i] = T[]
        end
        new(cells, D[1], D[2], D[3])
	end
end

KspaceGrid3{T}(cell::PBCCell, d::Real) where{T} =
        KspaceGrid3{T}(cell.inv * [d,d,d])

Base.size(g3::Grid3) = size(g3.cells)

Base.getindex(g3::Grid3, i::Int) = getindex(g3.cells, i)

Base.IndexStyle(::Type{<:Grid3}) = IndexLinear()

function Base.empty!(g3::Grid3)
    for cell in g3.cells
        empty!(cell)
    end
end

function findcell(g3::OrthoGrid3, (Ax,Ay,Az)::AbstractVector{<:Real})
    d = g3.d
    xn = floor(Int, (Ax - g3.Ox) / d) + 1
    yn = floor(Int, (Ay - g3.Oy) / d) + 1
    zn = floor(Int, (Az - g3.Oz) / d) + 1
    xn, yn, zn
end

function findcell(g3::KspaceGrid3, (Ax,Ay,Az)::AbstractVector{<:Real})
    xn = floor(Int, Ax / g3.Dx) + 1
    yn = floor(Int, Ay / g3.Dy) + 1
    zn = floor(Int, Az / g3.Dz) + 1
    xn, yn, zn
end

function Base.push!(g3::Grid3{T}, A::AbstractVector{<:Real}, v::T) where {T}
	push!(g3.cells[findcell(g3, A)...], v)
	g3
end

function Base.append!(g3::Grid3{T}, R::AbstractMatrix{<:Real},
        V::AbstractVector{T}) where {T}
    for i = 1:ncols(R)
        A = @view R[:,i]
        push!(g3, A, V[i])
    end
    g3
end

findnear(g3::Grid3{T}, A::AbstractVector{<:Real}) where {T} =
        findnear!(T[], g3, A)

function findnear!(dest::AbstractVector{T}, g3::OrthoGrid3{T},
        A::AbstractVector{<:Real}) where {T}
	empty!(dest)
	nx, ny, nz = size(g3.cells)
	xn, yn, zn = findcell(g3, A)
	for xi = xn-1:xn+1
		for yi = yn-1:yn+1
			for zi = zn-1:zn+1
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

function findnear!(dest::AbstractVector{T}, g3::KspaceGrid3{T},
        A::AbstractVector{<:Real}) where {T}
	empty!(dest)
	nx, ny, nz = size(g3.cells)
	xn, yn, zn = findcell(g3, A)
	for xi = xn-1:xn+1
		for yi = yn-1:yn+1
			for zi = zn-1:zn+1
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
