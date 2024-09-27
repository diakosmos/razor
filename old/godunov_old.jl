module Godunov
using LinearAlgebra
#using StaticArrays
"""
    For now, we assume a constant (fixed) viscosity "ν".
"""

ν = 3.01 / 3
s = 3.01 # relaxation time
q = 3*ν/s

"""
A planetary particulate disk is called a "ring", but this terminology is confusing.
"""
mutable struct ring
    N::Integer              # number of cells
    β::Float64              # power-law index for viscosity
    Δt::Float64             # time step
    X_::Array{Float64,1}    # positions (nodes)
    ΔX_::Array{Float64,1}   # cell sizes (cell)
    U_::Array{Float64,1}    # velocity (on nodes)
    Σ_::Array{Float64,1}    # surface density
    J_::Array{Float64,1}    # mass flux
    r₊::Array{Float64,1}    # 2-component eigenvector, right-travelling wave (relative to fluid) [NB: independent of position]
    r₋::Array{Float64,1}    # 2-component eigenvector, left-travelling wave (relative to fluid) [NB: independent of position]
    λ₊::Array{Float64,1}    # eigenvalues, right-travelling (relative to fluid) wave
    λ₋::Array{Float64,1}    # eigenvalues, left-travelling (relative to fluid) wave
    Λ₊R::Any                # diagonal matrix of eigenvalues, Right-travelling (absolute) right-travelling (relative)
    Λ₊L::Any                # diagonal matrix of eigenvalues, Left-travelling (absolute) right-travelling (relative)
    Λ₋R::Any                #
    Λ₋L::Any                #
    α⁺::Array{Float64,1}    #
    α⁻::Array{Float64,1}    #
    t::Float64              # time
end
function ring(N=5;β=0,xf=3.0,xi=0.01,vfield="standard")#10*xf/N)
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
        Λ₊R =  Diagonal(λ₊ .* (λ₊ .> 0))
        Λ₊L = -Diagonal(λ₊ .* (λ₊ .< 0))
        Λ₋R =  Diagonal(λ₋ .* (λ₋ .> 0))
        Λ₋L = -Diagonal(λ₋ .* (λ₋ .< 0))
        return (r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₋L)
    end
    """ in progress.... """
    function setU(U_,X_,vfield)
        if vfield=="standard"
            U_ = 1.0 ./ X_.^4
        elseif vfield=="still"
            U_ *= 0.0
#        elseif vfield =
        end
    end
    xtemp_ = collect(range(xi,xf,length=N+1))
    X_ = xtemp_[2:end]
    ΔX_ = diff(xtemp_)
    U_ = 1.0 ./ X_.^4
    #U_ *= 0.0;
    U_ .= -0.1
    Σ_ = ones(N); Σ_[1] = 0.0;
    Σ_[1:(2*N÷5)] .= 0.0; Σ_[(3*N÷5):end] .= 0.0
    J_ = zeros(N)
    (r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₋L) = setEs(U_)
    FUDGE = 0.1
    Δt = 0.5 / max( maximum(abs.(λ₊./ΔX_)), maximum(abs.(λ₋./ΔX_))) * FUDGE
    α⁺=zeros(N)
    α⁻=zeros(N)
    return ring(N,β,Δt,X_,ΔX_,U_,Σ_,J_,r₊,r₋,λ₊,λ₋,Λ₊R,Λ₊L,Λ₋R,Λ₋L,α⁺,α⁻,0.0)
end

"""
getR(): get right eigenvectors (two 2-component vectors), and eigenvalues (2*N).
"""


function getAlphas!(r::ring)
    N = r.N; Σ_ = r.Σ_; J_ = r.J_; α⁺ = r.α⁺; α⁻ = r.α⁻
    for i in 1:N
        if i < N
            ΔΣ = Σ_[i+1] - Σ_[i]
            ΔJ = J_[i+1] - J_[i]
        else
            ΔΣ = Σ_[1] - Σ_[i] # 1.0 - Σ_[i]
            ΔJ = J_[1] - J_[i] #0.0 - J_[i]
        end
        #print(0.5 * [1  1/√q] ⋅ [ΔΣ, ΔJ])
        α⁺[i] = 0.5 * [1  1/√q] ⋅ [ΔΣ, ΔJ]
        α⁻[i] = 0.5 * [1 -1/√q] ⋅ [ΔΣ, ΔJ]
    #=    α⁺[i] = [ √q , 1] ⋅ [ΔΣ, ΔJ] ./ √(1+q)
        α⁻[i] = [ √q ,-1] ⋅ [ΔΣ, ΔJ] ./ √(1+q) =#
    end
end#function

function step!(r::ring)
    getAlphas!(r)
    Δt=r.Δt; Σ_=r.Σ_; J_=r.J_; α⁺=r.α⁺; α⁻=r.α⁻; r₋=r.r₋; r₊=r.r₊
    Λ₊L=r.Λ₊L; Λ₊R=r.Λ₊R; Λ₋L=r.Λ₋L; Λ₋R=r.Λ₋R
    N = length(Σ_)
    """ LHS: """
    Σ_[1:end  ] += Δt * ((Λ₋L * α⁻ * r₋[1])                     )[1:end  ]
    Σ_[1:end  ] += Δt * (                     (Λ₊L * α⁺ * r₊[1]))[1:end  ]
    Σ_[2:  end] += Δt * (-(Λ₋R * α⁻ * r₋[1])                     )[1:end-1]
    Σ_[2:  end] += Δt * (                     -(Λ₊R * α⁺ * r₊[1]))[1:end-1]
    #Σ_ = Σ_.* (Σ_ .> 0)  # sanity - kills instability?
    J_[1:end  ] += Δt * ((Λ₋L * α⁻ * r₋[2])                     )[1:end  ]
    J_[1:end  ] += Δt * (                      (Λ₊L * α⁺ * r₊[2]))[1:end  ]
    J_[2:  end] += Δt * (-(Λ₋R * α⁻ * r₋[2])                     )[1:end-1]
    J_[2:  end] += Δt * (                     -(Λ₊R * α⁺ * r₊[2]))[1:end-1]
    #=
    #Threads.@threads for i in 1:N
    @simd for i in 1:N
        if i<N
            Σ_[i] += Δt * ((Λ₋L[i] * α⁻[i] * r₋[1])                     )
            Σ_[i] += Δt * (                     (Λ₊L[i] * α⁺[i] * r₊[1]))
        end
        if i>1
            Σ_[i] += Δt * (-(Λ₋R[i-1] * α⁻[i-1] * r₋[1])                     )
            Σ_[i] += Δt * (                     -(Λ₊R[i-1] * α⁺[i-1] * r₊[1]))
        end
        if i<N
            J_[i] += Δt * ((Λ₋L[i] * α⁻[i] * r₋[2])                     )
            J_[i] += Δt * (                     (Λ₊L[i] * α⁺[i] * r₊[2]))
        end
        if i>1
            J_[i] += Δt * (-(Λ₋R[i-1] * α⁻[i-1] * r₋[2])                     )
            J_[i] += Δt * (                     -(Λ₊R[i-1] * α⁺[i-1] * r₊[2]))
        end
    r.Σ_[i] = Σ_[i]
    r.J_[i] = J_[i]
    end
    =#
    """ RHS: """

    J_[1:end] -= Δt * J_[1:end] / s # forwards Euler

    r.t += Δt
end

end#module
