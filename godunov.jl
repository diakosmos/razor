module Godunov
using LinearAlgebra
#using StaticArrays
"""
    For now, we assume a constant (fixed) viscosity "ν".
"""

ν = 1.0
s = 1.0 # relaxation time
q = 3*ν/s

"""
A planetary particulate disk is called a "ring", but this terminology is confusing.
"""
mutable struct ring
    N::Integer # number of cells
    β::Float64 # power-law index for viscosity
    Δt::Float64
    X_::Array{Float64,1}
    ΔX_::Array{Float64,1}
    U_::Array{Float64,1}
    Σ_::Array{Float64,1}
    J_::Array{Float64,1}
    r₊::Array{Float64,1}
    r₋::Array{Float64,1}
    λ₊::Array{Float64,1}
    λ₋::Array{Float64,1}
    Λ₊R::Any
    Λ₊L::Any
    Λ₋R::Any
    Λ₋L::Any
    α⁺::Array{Float64,1}
    α⁻::Array{Float64,1}
end
function ring(N=5;β=0,xf=3.0,xi=0.01)#10*xf/N)
    xtemp_ = collect(range(xi,xf,length=N+1))
    X_ = xtemp_[2:end]
    ΔX_ = diff(xtemp_)
    U_ = 1.0 ./ X_.^4
    Σ_ = ones(N); Σ_[1] = 0.0
    J_ = zeros(N)
    (r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₊L) = setEs(U_)
    Δt = 0.5 / max( maximum(abs.(λ₊./ΔX_)), maximum(abs.(λ₋./ΔX_)))
    α⁺=zeros(N)
    α⁻=zeros(N)
    return ring(N,β,Δt,X_,ΔX_,U_,Σ_,J_,r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₊L,α⁺,α⁻)
end

"""
getR(): get right eigenvectors (two 2-component vectors), and eigenvalues (2*N).
"""
function setEs(U_::Array{Float64,1})
    N = length(U_)
    r₊ = [1.0, √q]
    r₋ = [1.0,-√q]
    λ₊ = zeros(N)
    λ₋ = zeros(N)
    for i in 1:N
        λ₊[i] = U_[i] + √q
        λ₋[i] = U_[i] - √q
    end
    Λ₊R = Diagonal(λ₊ .* (λ₊ .> 0))
    Λ₊L = Diagonal(λ₊ .* (λ₊ .< 0))
    Λ₋R = Diagonal(λ₋ .* (λ₋ .> 0))
    Λ₋L = Diagonal(λ₋ .* (λ₋ .< 0))
    return (r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₋L)
end

function getAlphas(r::ring)
    N = r.N; Σ_ = r.Σ_; J_ = r.J_; α⁺ = r.α⁺; α⁻ = r.α⁻
    for i in 1:N
        if i < N
            ΔΣ = Σ_[i+1] - Σ_[i]
            ΔJ = J_[i+1] - J_[i]
        else
            ΔΣ = 1 - Σ_[i]
            ΔJ = 0.0
        end
        print(0.5 * [1  1/√q] ⋅ [ΔΣ, ΔJ])
        α⁺[i] = 0.5 * [1  1/√q] ⋅ [ΔΣ, ΔJ]
        α⁻[i] = 0.5 * [1 -1/√q] ⋅ [ΔΣ, ΔJ]
    end
end#function

function step!(r::ring)
    getAlphas(r)
    Δt=r.Δt; Σ_=r.Σ_; α⁺=r.α⁺; α⁻=r.α⁻;
    Λ₊L=r.Λ₊L; Λ₊R=r.Λ₊R
    Σ_[1:end-1] += Δt * α⁻ .* λ₋
end

end#module
