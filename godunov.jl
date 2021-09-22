module Godunov
using LinearAlgebra
#using StaticArrays
"""
    For now, we assume a constant (fixed) viscosity "ν".
"""

TINY = 0.5 # how much to shink the timestep

α = 1.0
ν = 1.0
s = 1.0 # relaxation time
q = 3*ν/s
xv = (α/(3ν))^(1.0/3.0)
xh = (α/√q)^0.25
print("xv: $xv\n")
print("xh: $xh\n")
# Boundary conditions (Dirichlet):
Σ_lbc = 0.0
Σ_rbc = 1.0
Φ_lbc = 0.0
Φ_rbc = 0.0#-0.4

"""
A planetary particulate disk is called a "ring", but this terminology is confusing.
"""
mutable struct ring
    N::Integer              # number of cells
    Δt::Float64             # time step
    X_::Array{Float64,1}    # positions (nodes)
    ΔX_::Array{Float64,1}   # cell sizes (cell)
    U_::Array{Float64,1}    # velocity (on nodes)
    ∇U_::Array{Float64,1}   # divergence of velocity
    Σ_::Array{Float64,1}    # surface density
    S_::Array{Float64,1}    # integral of the above
    Φ_::Array{Float64,1}    # mass flux
    F_::Array{Float64,1}    # integral of the above
    r₊::Array{Float64,1}    # 2-component eigenvector, right-travelling wave (relative to fluid) [NB: independent of position]
    r₋::Array{Float64,1}    # 2-component eigenvector, left-travelling wave (relative to fluid) [NB: independent of position]
    λ₊::Array{Float64,1}    # eigenvalues, right-travelling (relative to fluid) wave
    λ₋::Array{Float64,1}    # eigenvalues, left-travelling (relative to fluid) wave
#    Λ₊R::Any                # diagonal matrix of eigenvalues, Right-travelling (absolute) right-travelling (relative)
#    Λ₊L::Any                # diagonal matrix of eigenvalues, Left-travelling (absolute) right-travelling (relative)
#    Λ₋R::Any                #
#    Λ₋L::Any                #
    α⁺::Array{Float64,1}    #
    α⁻::Array{Float64,1}    #
    t::Float64              # time
    i::Integer              # iteration counter
end
function ring(N=5;xi=(0.5*min(xv,xh)),xf=(3.0*max(xv,xh)))
    function setEs(U_::Array{Float64,1})
        n = length(U_)
        r₊ = [1.0, √q]
        r₋ = [1.0,-√q]
        λ₊ = zeros(n)
        λ₋ = zeros(n)
        for i in 1:n
            λ₊[i] = U_[i] + √q
            λ₋[i] = U_[i] - √q
        end
#=        Λ₊R =  Diagonal(λ₊ .* (λ₊ .> 0))
        Λ₊L = -Diagonal(λ₊ .* (λ₊ .< 0))
        Λ₋R =  Diagonal(λ₋ .* (λ₋ .> 0))
        Λ₋L = -Diagonal(λ₋ .* (λ₋ .< 0)) =#
        return (r₊,r₋,λ₊,λ₋)#,Λ₊R,Λ₊L,Λ₋R,Λ₋L)
    end

    x = 5.0; y = 1.0/x
    X_ = collect(range(xi^y,xf^y,length=N+1)).^x
    ΔX_ = diff(X_)
    U_ = α ./ X_.^4
    #U_ .= -1.5
    ∇U_ = diff(U_) ./ ΔX_

    #Σ_ = ones(N); Σ_[1] = 0.0; #Σ_[N÷2:end].=1.0
    Σ_ = (X_[2:end].-xh)./xf
    Σ_ .*= (Σ_.>0)
    Φ_ = Σ_ .* (U_[1:end-1]+U_[2:end])/2.0  .* 0.0

    S_ = Σ_ .* ΔX_
    F_ = Φ_ .* ΔX_

    (r₊,r₋,λ₊,λ₋) = setEs(U_)

    Δt = 0.5 / max( maximum(abs.(λ₊[1:end-1]./ΔX_)), maximum(abs.(λ₋[2:end]./ΔX_)), 1.0/s, maximum(abs.(∇U_))) * TINY
    print("Δt: $Δt\n")
    α⁺=zeros(N+1)
    α⁻=zeros(N+1)
    r = ring(N,Δt,X_,ΔX_,U_,∇U_,Σ_,S_,Φ_,F_,r₊,r₋,λ₊,λ₋,α⁺,α⁻,0.0,0)
    getAlphas!(r)
    return r
end

function getAlphas!(r::ring)
    N = r.N; Σ_ = r.Σ_; Φ_ = r.Φ_; α⁺ = r.α⁺; α⁻ = r.α⁻
    ΔΣ_ = cat(Σ_[1]-Σ_lbc, diff(Σ_), Σ_rbc-Σ_[end], dims=1)
    ΔΦ_ = cat(Φ_[1]-Φ_lbc, diff(Φ_), Φ_rbc-Φ_[end], dims=1)
    for i in 1:(N+1)
        α⁺[i] = 0.5 * [1  1/√q] ⋅ [ΔΣ_[i], ΔΦ_[i]]
        α⁻[i] = 0.5 * [1 -1/√q] ⋅ [ΔΣ_[i], ΔΦ_[i]]
    end
    r.α⁺ = α⁺
    r.α⁻ = α⁻
end#function

function step!(r::ring; Δt=r.Δt)
    N=r.N; Σ_=r.Σ_; ΔX_=r.ΔX_; S_=r.S_; Φ_=r.Φ_; F_=r.F_; α⁺=r.α⁺; α⁻=r.α⁻; r₋=r.r₋; r₊=r.r₊
    λ₊=r.λ₊; λ₋=r.λ₋
    getAlphas!(r)

    S_ = Σ_ .* ΔX_
    F_ = Φ_ .* ΔX_
    """ LHS: """
    # right-traveling.
    Σ_ -= r.α⁺[1:end-1] .* Δt .* r₊[1] .* λ₊[1:end-1] .* (λ₊[1:end-1].>0)
    Σ_ -= r.α⁺[1:end-1] .* Δt .* r₋[1] .* λ₋[1:end-1] .* (λ₋[1:end-1].>0)
    Φ_ -= r.α⁺[1:end-1] .* Δt .* r₊[2] .* λ₊[1:end-1] .* (λ₊[1:end-1].>0)
    Φ_ -= r.α⁺[1:end-1] .* Δt .* r₋[2] .* λ₋[1:end-1] .* (λ₋[1:end-1].>0)
    # left-traveling
    Σ_ -= r.α⁻[2:end]   .* Δt .* r₊[1] .* λ₊[2:end] .* (λ₊[2:end].<0)
    Σ_ -= r.α⁻[2:end]   .* Δt .* r₋[1] .* λ₋[2:end] .* (λ₋[2:end].<0)
    Φ_ -= r.α⁻[2:end]   .* Δt .* r₊[2] .* λ₊[2:end] .* (λ₊[2:end].<0)
    Φ_ -= r.α⁻[2:end]   .* Δt .* r₋[2] .* λ₋[2:end] .* (λ₋[2:end].<0)
    """ RHS: """
    # This part kinda sucks b/c it doesn't conserve very well:
    # Σ_ = S_ ./ ΔX_;  Φ_ = F_ ./ ΔX_
    Σ_[1] = Σ_lbc
    Σ_[end] = Σ_rbc
    Φ_[1] = Φ_lbc
    Φ_[end] = Φ_rbc

    #Σ_ += -Δt * r.∇U_ .* Σ_ # <- forward Euler

    Σ_ ./= (1 .+ Δt*r.∇U_/10)  # <- backward Euler  - conserves a bit better
    #Σ_ .*= (1 .- Δt*r.∇U_)
    #Φ_ ./= (1 .+ Δt*r.∇U_ .+ Δt/s)
    Φ_ ./= (1 .+ Δt/s)
#=
    for i in 1:N
        if r.X_[i]<xh
            Σ_[i] = Φ_[i] = 0
        end
    end
=#
    #


    #S_ *= 0
    #
    # massaging.
    # 1. make sign uniform
#    Σ_ = abs.(S_); Φ_ = -abs.(Φ_)
    # 2. make monotonic
    for i in 2:N
        Σ_[i] = max(Σ_[i-1],Σ_[i])
#        Φ_[i] = min(Φ_[i-1],Φ_[i])
    end
    Σ_ /= maximum(Σ_)
# =#
#    Σ_[7] = 0.0
    r.Σ_ = Σ_; r.Φ_ = Φ_
    S_ = Σ_ .* ΔX_; F_ = Φ_ .* ΔX_
    r.S_ = S_
#    r.Σ_ = S_ ./ ΔX_
    r.Φ_ = Φ_
    r.F_ = F_  #Φ_ ./ ΔX_

    r.t += Δt
    r.i += 1
    s1 = sum(S_); s2 = sum(F_)
    #print("Sums: $s1,  $s2 \n")
    return (sum(S_),sum(F_))
end

function tstep!(r::ring,Δt)
    M = floor(Δt / r.Δt)
    Δt_last = Δt - M*r.Δt
    @assert Δt_last ≥ 0
    for i in 1:M
        step!(r)
    end
    step!(r, Δt = Δt_last)
end

end#module
