module Sand_03
using LinearAlgebra, Random

#export ring, tstep!, step!

#using StaticArrays
"""
    Sand 03: Hyperbolic diffusion, with no advection (U=0)
"""

JIGGLE = 0.9 # how much to randomly tweak the node positions, to test algo; 0 ≤ JIGGLE < 1.0
q = 1.0 # hyperbolic wave speed
s = 3.0 # relaxation time
r₊ = [1.0, √q]
r₋ = [1.0,-√q]

mutable struct ring
    N::Integer              # number of cells
    Δt::Float64             # time step
    X_::Array{Float64,1}    # positions (nodes)
    ΔX_::Array{Float64,1}   # cell sizes (cell)
#    U_::Array{Float64,1}    # velocity (on nodes) XXX IGNORED
#    ∇U_::Array{Float64,1}   # divergence of velocity (on cells) XXX IGNORED
    Σ_::Array{Float64,1}    # surface density
    S_::Array{Float64,1}    # integral of the above
    J_::Array{Float64,1}    # Flux
    F_::Array{Float64,1}    # integral of flux
    α⁺::Array{Float64,1}    #
    α⁻::Array{Float64,1}    #
    t::Float64              # time
end
function ring(N=100; U = 1.0, xi=0.0, xf=10.0)
    xtemp_ = collect(range(xi,xf,length=N+1))
    xtemp_[2:end-1] += (Random.rand(N-1).-0.5) * (xf-xi)/N * JIGGLE
    @assert prod(diff(xtemp_).>0)  # checks to make sure xtemp_ monotonic
    X_ = xtemp_#[2:end]
    ΔX_ = diff(xtemp_)

    # The next several lines are superfluous since
#    U_ = U * ones(N+1)
#    """ XXX putting step-down in speed XXX """
#    U_[3*N÷4:end] .*= 0.5
#    U_[1:N÷4]     .*= 0.5
#    ∇U_ = diff(U_) ./ ΔX_
    # init w/top-hat density:
    Σ_ = ones(N); Σ_[1] = 0.0;
    Σ_[1:(2*N÷5)] .= 0.0; Σ_[(3*N÷5):end] .= 0.0
    S_ = Σ_ .* ΔX_
    J_ = zeros(N); F_ = J_ * 1.0
    #
    #Δt = 0.5 / max(maximum(U_[1:end-1]./ΔX_),-minimum(U_[2:end]./ΔX_))
    Δt = 0.5 / max( maximum(q./ΔX_), -minimum(-q./ΔX_) , s)
    #
    return ring(N,Δt,X_,ΔX_,Σ_,S_,J_,F_,zeros(N+1),zeros(N+1),0.0)
end

function getAlphas!(r::ring)
    N = r.N; Σ_ = r.Σ_; J_ = r.J_; α⁺ = r.α⁺; α⁻ = r.α⁻
    ΔΣ_ = cat(Σ_[1]-Σ_[end], diff(Σ_), Σ_[1]-Σ_[end], dims=1)
    ΔJ_ = cat(J_[1]-J_[end], diff(J_), J_[1]-J_[end], dims=1)
    for i in 1:(N+1) # Tried to write this as two one-liners instead of loop, but probs w/converting Vec{Array} to Array...
        α⁺[i] = 0.5 * [1  1/√q] ⋅ [ΔΣ_[i], ΔJ_[i]]
        α⁻[i] = 0.5 * [1 -1/√q] ⋅ [ΔΣ_[i], ΔJ_[i]]
    end
    r.α⁺ = α⁺
    r.α⁻ = α⁻
end#function

function step!(r::ring; Δt=r.Δt)
    Σ_=r.Σ_; N = r.N; ΔX_ = r.ΔX_; S_ = r.S_; J_=r.J_; F_ =r.F_
    #
    getAlphas!(r)
    #
    # get integral quantities
    S_ = Σ_ .* ΔX_
    F_ = J_ .* ΔX_
#    U⁺_ = U_[1:end-1] .* (U_[1:end-1].>0) # right-traveling
#    U⁻_ = U_[2:end  ] .* (U_[2:end  ].<0) # left-traveling
    # right-traveling. nb: √q is wave speed
    S_ -= r.α⁺[1:end-1] * Δt * r₊[1] * √q
    F_ -= r.α⁺[1:end-1] * Δt * r₊[2] * √q
    # left-traveling
    S_ += r.α⁻[2:end] * Δt * r₋[1] * √q
    F_ += r.α⁻[2:end] * Δt * r₋[2] * √q


    # This part kinda sucks b/c it doesn't conserve very well:
    Σ_ = S_ ./ ΔX_;  J_ = F_ ./ ΔX_
    #Σ_ += -Δt * r.∇U_ .* Σ_ # <- forward Euler
    #  Σ_ ./= (1 .+ Δt*r.∇U_)  # <- backward Euler  - conserves a bit better
    J_ ./= (1 + Δt / s)
    #
    S_ = Σ_ .* ΔX_; F_ = J_ .* ΔX_


    #
    r.S_ = S_
    r.Σ_ = S_ ./ ΔX_
    r.F_ = F_
    r.J_ = F_ ./ ΔX_

    r.t += Δt
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
