using SparseArrays
using PowerModels
using LinearAlgebra

abstract type AbstractPTDF end


struct LazyPTDF{TF} <: AbstractPTDF
    N::Int  # number of buses
    E::Int  # number of branches
    islack::Int  # Index of slack bus

    A::SparseMatrixCSC{Float64,Int}  # incidence matrix
    b::Vector{Float64}  # branch susceptances
    BA::SparseMatrixCSC{Float64,Int}  # B*A
    AtBA::SparseMatrixCSC{Float64,Int}  # AᵀBA

    F::TF   # Factorization of -(AᵀBA). Must be able to solve linear systems with F \ p
            # We use a factorization of -(AᵀBA) to support cholesky factorization when possible

    # TODO: cache
end

function LazyPTDF(data::OPFData; solver::Symbol=:ldlt)
    N, E, A, b, ref_idx = data.N, data.E, data.A, data.b, data.ref_bus

    B = Diagonal(b)
    BA = B * A
    S = AtBA = A' * BA

    S[ref_idx, :] .= 0.0
    S[:, ref_idx] .= 0.0
    S[ref_idx, ref_idx] = -1.0;  # to enable cholesky
    S = -S

    if solver == :ldlt
        F = ldlt(S)
    elseif solver == :cholesky
        # If Cholesky is not possible, default to LDLᵀ
        if maximum(b) < 0.0
            F = cholesky(S)
        else
            @warn "Some branches have positive susceptance, cannot use Cholesky; defaulting to LDLᵀ"
            F = ldlt(S)
        end
    else
        error("Invalid linear solver: only ldlt and cholesky are supported")
    end

    return LazyPTDF(N, E, ref_idx, A, b, BA, AtBA, F)
end

"""
    compute_flow_lazy!(pf, pg, Φ::LazyPTDF)

Compute power flow `pf = Φ*pg` lazyly, without forming the PTDF matrix.

Namely, `pf` is computed as `pf = BA * (F \\ pg)`, where `F` is a factorization
    of (-AᵀBA), e.g., a cholesky / LDLᵀ / LU factorization.
"""
function compute_flow!(pf, pg, Φ::LazyPTDF)
    θ = Φ.F \ pg
    θ[Φ.islack, :] .= 0  # slack voltage angle is zero
    mul!(pf, Φ.BA, θ, -one(eltype(pf)), zero(eltype(pf)))
    return pf
end

"""
    ptdf_row(Φ::LazyPTDF, e::Int)

Return the `e`-th row of (lazy) PTDF matrix `Φ`.
"""
function ptdf_row(Φ::LazyPTDF, e::Int)
    1 <= e <= Φ.E || error("Invalid row index: $e (must be between 1 and $E)")

    u = -Φ.BA[e, :]
    y = Φ.F \ u
    y[Φ.islack] = 0.0

    return y
end


struct FullPTDF{T} <: AbstractPTDF
    N::Int  # number of buses
    E::Int  # number of branches

    matrix::Matrix{T}  # PTDF matrix
end

function FullPTDF(data::OPFData)
    lazy = LazyPTDF(data)
    return FullPTDF(lazy)
end
function FullPTDF(lazy::LazyPTDF)
    T = eltype(lazy.F)
    PTDF = zeros(T, lazy.E, lazy.N)
    for e in 1:lazy.E
        PTDF[e, :] .= ptdf_row(lazy, e)
    end

    return FullPTDF(lazy.N, lazy.E, PTDF)
end

"""
    compute_flow_direct!(pf, pg, Φ::FullPTDF)

Compute power flow `pf = Φ*pg` given PTDF matrix `Φ` and nodal injections `pg`.
"""
function compute_flow!(pf, pg, Φ::FullPTDF)
    mul!(pf, Φ.matrix, pg)
    return pf
end

"""
    ptdf_row(Φ::FullPTDF, e::Int)

Return the `e`-th row of PTDF matrix `Φ`.
"""
function ptdf_row(Φ::FullPTDF, e::Int)
    return Φ.matrix[e, :]
end