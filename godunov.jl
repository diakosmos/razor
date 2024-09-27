module Godunov
"""
GodunovFull: "Full" means that we allow the eigenvectors to depend upon position,
b/c we are treating the "full" problem in which the viscosity ν (and therefore the diffusion D)
and the relaxation time s are not fixed, but depend upon density Σ, in the form of 
a power-law. This makes the problem more difficult, b/c the Godunov eigenvectors now 
depend upon Σ, which changes as a function of time and location.

Futhermore, we allow for the possibility that U is a function of time. The intent is to play
with either letting it be independent of time (as before), representing the time-averaged forcing
over a full relative orbit, or to try some method of incorporating the time-dependence of the forcing.
That's tricky, though, b/c the time between successive forcings depends on how far removed you are,
radially, from the moon.
"""

using Plots
using LinearAlgebra
using Random
Random.seed!(92842)

TINY = 0.5 # how much to shrink the timestep 1.0 is CFL (supposedly, anyway)

mutable struct ring
    N::Integer              # number of cells
    parmd::Dict{String,Float64} # fundamental physical parameter

    Ls::Dict{String,Float64}  # dictionary of fundamental lengthscales

    X●::Array{Float64,1}    # positions (nodes)  (N+1)
    ΔX_::Array{Float64,1}   # cell sizes (cell)  (N+2) (2 ghost cells)

    U●::Array{Float64,1}    # velocity (on nodes)
    # ∇U_::Array{Float64,1}   # divergence of velocity

    Σ_::Array{Float64,1}    # surface density
    S_::Array{Float64,1}    # integral of the above

    Φ_::Array{Float64,1}    # mass flux
    F_::Array{Float64,1}    # integral of the above

    r₊::Array{Float64,1}    # 2-component eigenvector, right-travelling wave (relative to fluid) [NB: dependent on time & position]
    r₋::Array{Float64,1}    # 2-component eigenvector, left-travelling wave (relative to fluid) [NB: dependent on time & position]

    λ₊●::Array{Float64,1}    # eigenvalues, right-travelling (relative to fluid) wave
    λ₋●::Array{Float64,1}    # eigenvalues, left-travelling (relative to fluid) wave

    α⁺●::Array{Float64,1}    # right-traveling amplitude in eigenvector decomposition
    α⁻●::Array{Float64,1}    # left

    Δt::Float64             # time step
    t::Float64              # time counter
    i::Integer              # iteration counter
end
function ring(N=20, params=(1.0,1.0/3,1.0); xf = 3.0, xexp=1)
    # Get key physical parameters, and lengthscales to non-dimensionalize X-axis.
    α = params[1]; ν = params[2]; s = params[3] # velocity amplitude, kinematic viscosity, relaxation time
    q = 3*ν/s; D = 3*ν;    X = α^2 / (3ν)^5 / s^3
    parmd = Dict("α"=>α, "ν"=>ν, "s"=>s, "q"=>q, "D"=>D, "X"=>X)
    #
    Lc = (α/√q)^0.25        # position of critical point for hyperbolic transport (could call it "L0" I suppose)
    Lv = (α/(3ν))^(1.0/3.0) # viscous "critical" point
    L2 = (α*s)^0.2          # another lengthscale constructed from α, ν, s
    L3 = √(ν*s)             # another lengthscale      "        "    "   "
 
    Ls = Dict("Lcrit"=>Lc, "Lvisc"=>Lv, "L2"=>L2, "L3"=>L3)
    a = sort(collect(Ls),by=x->x[2]) # this is a bit ugly... ugh
    for e in a
        print("$e\n")
    end
    #
    # Set X positions (nodes):
    x = xexp; y=1.0/x # exponents to space-out nodes by power law
    xi = Lc
    X● = collect(range(xi^y,xf^y,length=N+1)).^x
    ΔX_ = diff(X●)
    ΔX_ = cat(ΔX_[1],ΔX_,ΔX_[end],dims=1) # add ghost cells on either end
    #
    # Set speed and divergence of same
    X2● = 2 * X●[end] .- X● 
    U● = α ./ X●.^4 - α ./ X2●.^4   
    #∇U_ = diff(U_) ./ ΔX_
    #
    # Initialize mass density by making it uniform but zero for X < L_c
    Σ_ = ones(N+2) # w/ghost cell on either end
    for i in 1:N+1
        if X●[i] < Ls["Lcrit"]
            Σ_[i] = 0.0
        else 
            Σ_[i] = 2*Random.rand() 
        end
    end
    # BCs:
    Σ_[1] = 0.0; Σ_[2] = 1.0; Σ_[end]=Σ_[end-1] = 1.0
    S_ = Σ_ .* ΔX_

    # Initialize flux to suface density times velocity (to reduce initial transients)
    Φ_ = zeros(N+2)
    Φ_[2:end-1] = Σ_[2:end-1] .* U●[2:end] # just a no-too-horrible initial guess
    Φ_[2] = -Σ_[2] * √q ;  Φ_[end] = -Φ_[end-1]
    F_ = Φ_ .* ΔX_ # it is odd to think of flux as a conserved quantity, but it can be useful at times...
    #
    # eigenvectors, eigenvalues
    #r₊● = [zeros(2) for _ in 1:N+1]
    r₊ = [1.0, √q]
    r₋ = [1.0,-√q]
    λ₊● = zeros(N+1)
    λ₋● = zeros(N+1)
    for i in 1:N+1
        λ₊●[i] = U●[i] + √q
        λ₋●[i] = U●[i] - √q
    end
    #
    # largest suggested timestep:
    Δt = 0.5 / max( maximum(abs.(λ₊●[1:end-1]./ΔX_[2:end-1])),
        maximum(abs.(λ₋●[2:end]./ΔX_[2:end-1])), 1.0/s) * TINY
    print("Δt: $Δt\n")
    #
    # amplitudes of eigenvector decomposition of ΔΣ, ΔΦ fields.
    α⁺●=zeros(N+1) 
    α⁻●=zeros(N+1)
    r = ring(N,parmd,Ls,X●,ΔX_,U●,Σ_,S_,Φ_,F_,r₊,r₋,λ₊●,λ₋●,
        α⁺●,α⁻●,Δt,0.0,0)
    getAlphas!(r)
    #
    # done.
    return r
end

function getAlphas!(r::ring)
    N=r.N; Σ_ = r.Σ_; Φ_ = r.Φ_; α⁺● = r.α⁺●; α⁻● = r.α⁻●;
    q = r.parmd["q"]
    #=ΔΣ● = cat(Σ_[1]-0.0, diff(Σ_),        0.0, dims=1) # reflective bc on rt
    ΔΦ● = cat(Φ_[1]-0.0, diff(Φ_), -Φ_[end], dims=1) # reflective bc w/sign flip
    for k in 1:(N+1)
        Δα⁺●[k] = 0.5 * [1  1/√q] ⋅ [ΔΣ●[k], ΔΦ●[k]]
        Δα⁻●[k] = 0.5 * [1 -1/√q] ⋅ [ΔΣ●[k], ΔΦ●[k]]
    end =#
    @simd for j in 1:N+1
        α⁺●[j] = 0.5 * [1  1/√q] ⋅ [Σ_[j+1]-Σ_[j],Φ_[j+1]-Φ_[j]]
        α⁻●[j] = 0.5 * [1 -1/√q] ⋅ [Σ_[j+1]-Σ_[j],Φ_[j+1]-Φ_[j]]
    end
    # BCs:
    #α⁺●[1] = 0.0;  #α⁻●[1] = 1.0
    α⁺●[N+1] = 0.5 * [1  1/√q] ⋅ [0, -2*Φ_[N+1]] # Sigma is even, Phi odd around bondary
    α⁻●[N+1] = 0.5 * [1 -1/√q] ⋅ [0, -2*Φ_[N+1]]

    r.α⁺● = α⁺●
    r.α⁻● = α⁻●
end#function

#=
"""
reconstruct!(): builds up Σ_ and Φ_ from α decomposision
"""
function reconstruct!(r::ring)#,α,β)
    N=r.N; Σ_ = r.Σ_; Φ_ = r.Φ_; α⁺● = r.α⁺●; α⁻● = r.α⁻●; S_ = r.S_; F_ = r.F_
    @simd for k in 2:N+1
        Σ_[k] = α⁺_[k] * r.r₊[1] + α⁻_[k] * r.r₋[1]
        Φ_[k] = α⁺_[k] * r.r₊[2] + α⁻_[k] * r.r₋[2]
    end
    Σ_[1] = Φ_[1] = 0.0
    Σ_[N+2] =  Σ_[N+1]
    Φ_[N+2] = -Φ_[N+1]
    if CROP
        for j in 2:N+1
            if r.X●[j-1] < r.Ls["Lcrit"]*EPSp1
                Σ_[j] = 0.0
                Φ_[j] = 0.0
            end
        end #
    end
    S_ = Σ_ .* r.ΔX_
    F_ = Φ_ .* r.ΔX_
    r.Σ_ = Σ_; r.Φ_ = Φ_; r.S_ = S_; r.F_ = F_
end =#

function step!(r::ring; Δt=r.Δt)
    N=r.N; Σ_=r.Σ_; ΔX_=r.ΔX_; S_=r.S_; Φ_=r.Φ_; F_=r.F_;
    λ₊●=r.λ₊●; λ₋●=r.λ₋●; α⁺● = r.α⁺●; α⁻● = r.α⁻●;
    X● = r.X●

    s = r.parmd["s"]

    getAlphas!(r)

    """ LHS: """
    #Σ_[j] -= α⁺[j] * Δt * r₋[1] * λ₋[j] .* (λ₋[j].>0) # <--- Example of another way to write the lines below.
    S_ = Σ_ .* ΔX_;  F_ = Φ_ .* ΔX_
    for i in 2:N+1
        S_[i] -= α⁻●[i-1] * λ₋●[i-1] * Δt * r.r₋[1] * (λ₋●[i-1]>0)
        S_[i] -= α⁺●[i-1] * λ₊●[i-1] * Δt * r.r₋[1] * (λ₊●[i-1]>0)
        S_[i] -= α⁻●[ i ] * λ₋●[ i ] * Δt * r.r₋[1] * (λ₋●[ i ]<0)
        S_[i] -= α⁺●[ i ] * λ₊●[ i ] * Δt * r.r₊[1] * (λ₊●[ i ]<0)
        F_[i] -= α⁻●[i-1] * λ₋●[i-1] * Δt * r.r₋[2] * (λ₋●[i-1]>0)
        F_[i] -= α⁺●[i-1] * λ₊●[i-1] * Δt * r.r₋[2] * (λ₊●[i-1]>0)
        F_[i] -= α⁻●[ i ] * λ₋●[ i ] * Δt * r.r₋[2] * (λ₋●[ i ]<0)
        F_[i] -= α⁺●[ i ] * λ₊●[ i ] * Δt * r.r₊[2] * (λ₊●[ i ]<0)
    end
    Σ_ = S_ ./ ΔX_;   Φ_ = F_ ./ ΔX_ 
    
    #for j in 2:N+1
    #    α⁺_[j] += (Δt/ΔX_[j]) * (-(λ₊●[ j ]>0 ? α⁺_[ j ] : α⁺_[j+1]) * λ₊●[ j ] +
    #                              (λ₊●[j-1]>0 ? α⁺_[j-1] : α⁺_[ j ]) * λ₊●[j-1] )
    #    α⁻_[j] += (Δt/ΔX_[j]) * (-(λ₋●[ j ]>0 ? α⁻_[ j ] : α⁻_[j+1]) * λ₋●[ j ] +
    #                              (λ₋●[j-1]>0 ? α⁻_[j-1] : α⁻_[ j ]) * λ₋●[j-1] )
    #end

    #r.α⁺_ = α⁺_
    #r.α⁻_ = α⁻_
    #reconstruct!(r)#,α⁺_,α⁻_)

    """ RHS: """
    ∇U_ = diff(r.U●) ./ ΔX_[2:end-1]; ∇U_ = cat(∇U_[1],∇U_,∇U_[end],dims=1)
    @assert maximum(abs.(∇U_))*Δt < 1
    Σ_ .*= exp.(-Δt*∇U_)
    Φ_ .*= exp.(-Δt*(∇U_.+1/s))

    #= or we could use analytic solution LOL
    X_ = X●[1:end-1] + X●[2:end]
    factor_ = 4 * r.parmd["α"] ./ X_.^3
    Σ_[2:end-1] .*= (1 .+ Δt .* factor_)
    Φ_[2:end-1] .*= (1 .+ Δt .* factor_) =#

    # Finish up, renormalize
    Σ_ .*= (Σ_.≥0)
    #Σ_[1] = 0.0; Φ_[1] = 0.0
    #Σ_[1]=1.0; Φ_[1] = -1.0/√r.parmd["q"]
    #
    # renormalize to make Σ_[end] = 1
    Σ_[end] = Σ_[end-1]; Φ_[end] = -Φ_[end-1]
    buildup = Σ_[end] 
    Σ_ ./= buildup; Φ_ ./= buildup
    S_ = Σ_ .* ΔX_; F_ = Φ_ .* ΔX_

    #
    # and increment the counters and return
    r.t += Δt
    r.i += 1
    #
    #mfig = plot(0.5*(r.X●[1:end-1]+r.X●[2:end]),r.Σ_[2:end-1])
    #display(mfig)
    r.Σ_ = Σ_; r.S_ = S_; r.Φ_ = Φ_; r.F_ = F_    
end


function tstep!(r::ring,Δt)
    M = floor(Δt / r.Δt)
    Δt_last = Δt - M*r.Δt
    #@assert Δt_last ≥ 0
    for i in 1:M
        step!(r)
    end
    step!(r, Δt = Δt_last)
end

end#module
