module GodunovFull
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

CFL_REL = 0.1 # time step relative to CFL condition; 0 < CFL_REL ≤ 1.0
TINY = 0.0 #1.0e-6 # stabilize when q = 0
Σmin = 0.01
CROP = true 
EPSp1 = 1.1
XEXP = 2

using LinearAlgebra

mutable struct ring
    N::Integer              # number of cells (plus 2 ghost cells)
    parmd::Dict{String,Float64} # fundamental physical parameters

    Ls::Dict{String,Float64}  # dictionary of fundamental lengthscales

    ν::Any # function 
    s::Any 
    q::Any 

    X●::Array{Float64,1}    # positions (nodes)
    ΔX_::Array{Float64,1}   # cell sizes (cell)

    U●::Array{Float64,1}    # velocity (on nodes)
    ΔU_::Array{Float64,1}   # difference of velocity

    Σ_::Array{Float64,1}    # surface density
    S_::Array{Float64,1}    # integral of the above

    Φ_::Array{Float64,1}    # mass flux
    F_::Array{Float64,1}    # integral of the above

    r₊●::Array{Float64,2}    # 2-component eigenvector, right-travelling wave (relative to fluid) [NB: dependent on time & position]
    r₋●::Array{Float64,2}    # 2-component eigenvector, left-travelling wave (relative to fluid) [NB: dependent on time & position]

    λ₊●::Array{Float64,1}    # eigenvalues, right-travelling (relative to fluid) wave
    λ₋●::Array{Float64,1}    # eigenvalues, left-travelling (relative to fluid) wave

    α⁺●::Array{Float64,1}    # right-traveling amplitude in eigenvector decomposition
    α⁻●::Array{Float64,1}    # left

    #β⁺_::Array{Float64,1}
    #β⁻_::Array{Float64,1}

    Δt::Float64             # time step
    t::Float64              # time counter
    i::Integer              # iteration counter
end

Keeler = Dict(
    "α" => 1.009e11, # m^5/s
    "ν" => 0.0250,   # m^2/s
    "s" => 2π * 8196.7,   # 2π/Ω
    "hΔx" => 18.5e3   # half gap
)

function ring(N=224, params=(Keeler["α"],Keeler["ν"],Keeler["s"]/100,1.0,0,0); xir = 0.95, xfr = 3.0) # = 0.95 * KeelerHalfGap, xf = 3.0*KeelerHalfGap)
    # Get key physical parameters, and lengthscales to non-dimensionalize X-axis.
    α = params[1]; ν0 = params[2]; s0 = params[3] # velocity amplitude, kinematic viscosity, relaxation time
    Σ0 = params[4]; β = params[5]; γ = params[6]
    q0 = 3*ν0/s0; D0 = 3*ν0
    parmd = Dict("α"=>α, "ν0"=>ν0, "s0"=>s0, "q0"=>q0, "D0"=>D0, "Σ0" => Σ0, "β" => β, "γ" => γ)
    #
    Lc = (α/√q0)^0.25        # position of critical point for hyperbolic transport (could call it "L0" I suppose)
    Lv = (α/(3ν0))^(1.0/3.0) # viscous "critical" point
    L2 = (α*s0)^0.2          # another lengthscale constructed from α, ν, s
    L3 = √(ν0*s0)             # another lengthscale      "        "    "   "
    Ls = Dict("Lcrit"=>Lc, "Lvisc"=>Lv, "L2"=>L2, "L3"=>L3)
    a = sort(collect(Ls),by=x->x[2]) # this is a bit ugly... ugh
    for e in a
        print("$e\n")
    end

    ν(Σ) = Σ>0 ? ν0 * (max(Σ,Σmin)/Σ0)^β : 0.0
    s(Σ) = Σ>0 ? s0 * (max(Σ,Σmin)/Σ0)^(-γ) : Inf
    q(Σ) = Σ>0 ? 3 * ν(Σ) / s(Σ) : 0.0
    #
    # Set X positions (nodes):
    x = XEXP; y=1.0/x # exponents to space-out nodes by power law
    xi = xir *Lv; xf = xfr*Lv
    X● = collect(range(xi^y,xf^y,length=N+1)).^x
    ΔX_ = cat(diff(X●)[1],diff(X●),diff(X●)[end],dims=1) # includes ghost cells
    #
    # Set speed and divergence of same
    X2● = 2 * X●[end] .- X●
    U● = α ./ X●.^4 - α ./ X2●.^4  # make U=0 on right boundary - symmetry
    ΔU_ = diff(U●) # ./ ΔX_
    #
    # Initialize mass density
    Σ_ = Σ0 * ones(N+2) # w/ghost cell on either end
    Σ_[1:N÷2] .= 0.0
    X_ = cat(X●[1],X●;dims=1) # just need this for initing Σ - nothing else
    #Σ_[X_.<Lv] .= 0.0

    S_ = Σ_ .* ΔX_
    #
    # Initialize flux to suface density times velocity (to reduce initial transients)
    Φ_ = zeros(N+2)
    Φ_[2:end-1] = Σ_[2:end-1] .* U●[2:end] # just a no-too-horrible initial guess
    F_ = Φ_ .* ΔX_ # it is odd to think of flux as a conserved quantity, but it can be useful at times...
    #
    # eigenvectors, eigenvalues
    r₊● = zeros(2,N+1) #[1.0, √q]
    r₋● = zeros(2,N+1) #[1.0,-√q]
    λ₊● = zeros(N+1)
    λ₋● = zeros(N+1)
    @simd for i in 1:N+1
        q₊ = q(Σ_[i+1]); q₋ = q(Σ_[i]) 
        qq = (q₊ + q₋)/2
        r₊●[1,i] = 1.0 ; r₋●[1,i] = 1.0
        r₊●[2,i] = √qq  ; r₋●[2,i] = -√qq
        λ₊●[i] = U●[i] + √qq
        λ₋●[i] = U●[i] - √qq
    end
    #
    # largest suggested timestep:
    Δt = 0.5 / max( maximum(abs.(λ₊●[1:end-1]./ΔX_[2:end-1])),
        maximum(abs.(λ₋●[2:end]./ΔX_[2:end-1])), maximum(1.0./s.(Σ_))) * CFL_REL
    print("Δt: $Δt\n")
    #
    # amplitudes of eigenvector decomposition of ΔΣ, ΔΦ fields.
    α⁺●=zeros(N+1) 
    α⁻●=zeros(N+1)
    r = ring(N,parmd,Ls,ν,s,q,X●,ΔX_,U●,ΔU_,Σ_,S_,Φ_,F_,r₊●,r₋●,λ₊●,λ₋●,
        α⁺●,α⁻●,Δt,0.0,1)
    getAlphas!(r)
    #
    # done.
    return r
end

function getAlphas!(r::ring)
    N=r.N; Σ_ = r.Σ_; Φ_ = r.Φ_; α⁺● = r.α⁺●; α⁻● = r.α⁻●;
    q = r.q
    ΔΣ● = diff(Σ_)
    ΔΦ● = diff(Φ_)
    @simd for j in 1:N+1
        qj = 0.5 * (q(Σ_[j+1])+q(Σ_[j]))
        qj = max(TINY,qj) # if qj = 0 it causes problems below:
        α⁺●[j] = 0.5 * [1  1/√qj] ⋅ [ΔΣ●[j],ΔΦ●[j]]
        α⁻●[j] = 0.5 * [1 -1/√qj] ⋅ [ΔΣ●[j],ΔΦ●[j]]
    end
    #α⁺●[1] = α⁻●[1] = 0.0
    #if CROP
    #    for j in 2:N+1
    #        if r.X●[j-1] < r.Ls["Lcrit"]*EPSp1
    #            α⁻●[j] = 0.0
    #        end
    #    end 
    #end
    r.α⁺● = α⁺●
    r.α⁻● = α⁻●
end

#="""
reconstruct!(): builds up Σ_ and Φ_ from α decomposision
"""
function reconstruct!(r::ring)#,α,β)
    N=r.N; Σ_ = r.Σ_; Φ_ = r.Φ_; α⁺_ = r.α⁺_; α⁻_ = r.α⁻_; S_ = r.S_; F_ = r.F_
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
    λ₊●=r.λ₊●; λ₋●=r.λ₋●; α₊● = r.α⁺●; α₋● = r.α⁻●; U● = r.U●
    X● = r.X●; r₋●=r.r₋●; r₊●=r.r₊●; ΔU_ = r.ΔU_; q = r.q

    s = r.s; ν = r.ν

    if r.i % 10 == 1 
    @simd for i in 1:N+1
        q₊ = q(Σ_[i+1]); q₋ = q(Σ_[i]) 
        qq = (q₊ + q₋)/2
        r₊●[1,i] = 1.0 ; r₋●[1,i] = 1.0
        r₊●[2,i] = √qq  ; r₋●[2,i] = -√qq
        λ₊●[i] = U●[i] + √qq
        λ₋●[i] = U●[i] - √qq
    end
    end
    r.r₊● = r₊●; r.r₋● = r₋●; r.λ₊● = λ₊●; r.λ₋● = λ₋●

    getAlphas!(r)

    #Δt = 0.5 / max( maximum(abs.(λ₊●[1:end-1]./ΔX_[2:end-1])),
    #    maximum(abs.(λ₋●[2:end]./ΔX_[2:end-1])), maximum(1.0./s.(Σ_))) * CFL_REL

    Δt = 0.5 / max( maximum(abs.(λ₊●[1:end-1]./ΔX_[2:end-1])),
         maximum(1.0./s.(Σ_))) * CFL_REL

    """ LHS: """
    @simd for j in 2:N+1
        S_[j] += Δt * ( (λ₋●[j-1]>0 ? -α₋●[j-1] * λ₋●[j-1] * r₋●[1,j-1] : 0.0) +
                        (λ₊●[j-1]>0 ? -α₊●[j-1] * λ₊●[j-1] * r₊●[1,j-1] : 0.0) +
                        (λ₋●[ j ]<0 ? α₋●[ j ] * λ₋●[ j ] * r₋●[1, j ] : 0.0) +
                        (λ₊●[ j ]<0 ? α₊●[ j ] * λ₊●[ j ] * r₊●[1, j ] : 0.0) )
        F_[j] += Δt * ( (λ₋●[j-1]>0 ? -α₋●[j-1] * λ₋●[j-1] * r₋●[2,j-1] : 0.0) +
                        (λ₊●[j-1]>0 ? -α₊●[j-1] * λ₊●[j-1] * r₊●[2,j-1] : 0.0) +
                        (λ₋●[ j ]<0 ? α₋●[ j ] * λ₋●[ j ] * r₋●[2, j ] : 0.0) +
                        (λ₊●[ j ]<0 ? α₊●[ j ] * λ₊●[ j ] * r₊●[2, j ] : 0.0) )
    end

    """ RHS: """
    S_[2:end-1] .+= -Δt .* ΔU_ .* Σ_[2:end-1]
    F_[2:end-1] .+= -Δt .* ΔU_ .* Φ_[2:end-1]
    s_ = s.(Σ_)
    Δν_ = cat(ν(Σ_[2])-ν(Σ_[1]),  ν.(Σ_[2:end]) .- ν.(Σ_[1:end-1]) ; dims=1 )
    F_ .-= Δt .* (3 .* Δν_ ./ s_) .* Σ_
    F_ ./= (1 .+ Δt./s_)
    # Finish up, renormalize
    Σ_ = S_ ./ ΔX_;  Φ_ = F_ ./ ΔX_
    Σ_ .*= (Σ_.≥0)
    if CROP
        Σ_[1] = 0.0; Φ_[1] = 0.0
        for j in 2:N+1
            if λ₋●[j-1] > 0  #r.X●[j-1] < r.Ls["Lcrit"]*EPSp1
                Σ_[j] = 0.0
                Φ_[j] = 0.0
            end
        end # =#
    end
    for j in 1:N+2
        Σ_[j] = max(Σmin,Σ_[j])
    end
    #Σ_[end-1:end] .= r.parmd["Σ0"]
    #Σ_[2:end-1] = 0.94 * Σ_[2:end-1] .+ 0.03 * (Σ_[1:end-2] + Σ_[3:end])
    S_ = Σ_ .* ΔX_; F_ = Φ_ .* ΔX_
    S_[1] = 0.0; Σ_[1] = 0.0; F_[1] = 0.0; Φ_[1] = 0.0
    F_[end] = 0.0; Φ_[end] = 0.0; F_[end-1] = 0.0; Φ_[end-1]=0.0
    r.S_ = S_; r.F_ = F_; r.Σ_ = Σ_; r.Φ_ = Φ_
    #
    # and increment the counters and return
    r.i += 1
    r.t += Δt
    return Δt
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
