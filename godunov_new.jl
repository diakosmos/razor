module Godunov
using LinearAlgebra

TINY = 0.5 # how much to shrink the timestep 1.0 is CFL (supposedly, anyway)

mutable struct ring
    N::Integer              # number of cells
    params::Tuple{Float64,Float64,Float64} # fundamental physical parameter

    Ls::Dict{String,Float64}  # dictionary of fundamental lengthscales

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

    α⁺::Array{Float64,1}    # right-traveling amplitude in eigenvector decomposition
    α⁻::Array{Float64,1}    # left

    Δt::Float64             # time step
    t::Float64              # time counter
    i::Integer              # iteration counter
end
function ring(N=20, params=(1.0,1.0,1.0); xir = 0.5, xfr = 3.0)
    # Get key physical parameters, and lengthscales to non-dimensionalize X-axis.
    α = params[1]; ν = params[2]; s = params[3] # velocity amplitude, kinematic viscosity, relaxation time
    Lc = (α/√q)^0.25        # position of critical point for hyperbolic transport (could call it "L0" I suppose)
    Lv = (α/(3ν))^(1.0/3.0) # viscous "critical" point
    L2 = (α*s)^0.2          # another lengthscale constructed from α, ν, s
    L3 = √(ν*s)             # another lengthscale      "        "    "   "
    Ls = Dict("Lcrit"=>Lc, "Lvisc"=>Lv, "L2"=>L2, "L3"=>L3)
    a = sort(collect(Ls),by=x->x[2]) # this is a bit ugly... ugh
    for e in a:
        print("$e[1]:  $e[2]\n")
    end
    #
    # Set X positions (nodes):
    x = 5.0; y=1.0/x # exponents to space-out nodes by power law
    X_ = collect(range(xi^y,xf^y,length=N+1)).^x
    ΔX_ = diff(X_)
    #
    # Set speed and divergence of same
    U_ = α ./ X_.^4;    U_ .-= U_[end]  # subtract off tiny bit, uniformly, to make U=0 on right boundary
    ∇U_ = diff(U_) ./ ΔX_
    #
    # Initialize mass density by making it uniform, with total mass = 1.0
    Σ_ = ones(N)
    S_ = Σ_ .* ΔX_
    S_ ./= sum(S_)
    Σ_ = S_ ./ ΔX_
    #
    # Initialize flux to suface density times velocity (to reduce initial transients)
    Φ_ = Σ_ .* U_[2:end] # U_[1:end-1] would have been fine too - just an initial guess
    F_ = Φ_ .* ΔX_ # it is odd to think of flux as a conserved quantity, but it can be useful at times...
    #
    # eigenvectors, eigenvalues
    r₊ = [1.0, √q]
    r₋ = [1.0,-√q]
    λ₊ = zeros(N+1)
    λ₋ = zeros(N+1)
    for i in 1:n
        λ₊[i] = U_[i] + √q
        λ₋[i] = U_[i] - √q
    end
    #
    # largest suggested timestep:
    Δt = 0.5 / max( maximum(abs.(λ₊[1:end-1]./ΔX_)), maximum(abs.(λ₋[2:end]./ΔX_)), 1.0/s, maximum(abs.(∇U_))) * TINY
    print("Δt: $Δt\n")
    #
    # amplitudes of eigenvector decomposition of Σ, Φ fields.
    α⁺=zeros(N+1)
    α⁻=zeros(N+1)
    r = ring(N,Ls,X_,ΔX_,U_,∇U_,Σ_,S_,Φ_,F_,r₊,r₋,λ₊,λ₋,α⁺,α⁻,Δt,0.0,0)
    getAlphas!(r)
    #
    # done.
    return r
end
