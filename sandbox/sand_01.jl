module Sand_01
using LinearAlgebra, Random
#using StaticArrays
"""
    Sand 01: simple uniform advection of a scalar
"""

mutable struct ring
    N::Integer              # number of cells
    Δt::Float64             # time step
    X_::Array{Float64,1}    # positions (nodes)
    ΔX_::Array{Float64,1}   # cell sizes (cell)
    U_::Array{Float64,1}    # velocity (on nodes)
    Σ_::Array{Float64,1}    # surface density
    S_::Array{Float64,1}    # integral of the above
    t::Float64              # time
end
function ring(N=100; U = 1.0, xi=0.0, xf=10.0)
    xtemp_ = collect(range(xi,xf,length=N+1))
    xtemp_[2:end-1] += (Random.rand(N-1).-0.5) * (xf-xi)/N * 0.9 # 0.8 is somewhat arbitrary
    @assert prod(diff(xtemp_).>0)  # checks to make sure xtemp_ monotonic
    X_ = xtemp_#[2:end]
    ΔX_ = diff(xtemp_)
    U_ = U * ones(N+1)
    """ U_[3*N÷4:end] .*= 0.5 # <--- this does NOT work properly; next step is to fix this. """
    # init w/top-hat density:
    Σ_ = ones(N); Σ_[1] = 0.0;
    Σ_[1:(2*N÷5)] .= 0.0; Σ_[(3*N÷5):end] .= 0.0
    S_ = Σ_ .* ΔX_
    #
    Δt = 0.5 / max(maximum(U_[1:end-1]./ΔX_),-minimum(U_[2:end]./ΔX_))
    #
    return ring(N,Δt,X_,ΔX_,U_,Σ_,S_,0.0)
end

function step!(r::ring; Δt=r.Δt)
    Σ_=r.Σ_; N = r.N; ΔX_ = r.ΔX_; U_ = r.U_; S_ = r.S_
    #
    ΔΣ_ = cat(Σ_[1]-Σ_[end], diff(Σ_), Σ_[1]-Σ_[end], dims=1)
    S_ = Σ_ .* ΔX_
    U⁺_ = U_[1:end-1] .* (U_[1:end-1].>0) # right-traveling
    U⁻_ = U_[2:end  ] .* (U_[2:end  ].<0) # left-traveling
    S_ -= Δt * U⁺_ .* ΔΣ_[1:end-1]
    S_ -= Δt * U⁻_ .* ΔΣ_[2:end]
    #
    r.S_ = S_
    r.Σ_ = S_ ./ ΔX_
    r.t += Δt
#    print(sum(S_),"\n")
    return sum(S_)
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
