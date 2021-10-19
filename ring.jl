module rings
using LinearAlgebra, DSP
include("cheb.jl")
using .cheb
export ring

"""
UNITS: cgs

"""

ASSERT = true # for debugging purposes: include "@assert" statements

# Consts
G  = 6.67430e-8 # Gravitational constant, cgs
M♄ = 5.6834e29  # mass of Saturn, in grams
A₁ = 6.7187     # Goldreich-Tremaine constant
Mm = 1.0e8      # Megameter. Typical orbit is 100Mm.

# Some fiducial satellites, for reference:
Pan  = (4.95e19, 1.33584e10) # (mass, orbital radius) of Pan (Encke gap)
Daph = (8.40e16, 1.36504e10) # (mass, orbital radius) of Daphnis (Keeler gap)

"""
"aring" is an abstract, idealized 1D (in r, i.e. in x) ring, with periodic b.c.
in the approximation that Δr ≪ r, and with a single mass forcing a gap.

- velocity is normalized to the diffusive wave speed √(κ₀/s₀)
- position is normalized to the critical lengthscale, Lcrit
- density Σ is normalized to the fiducial density Σ₀:  σ ≡ Σ/Σ₀
- the periodic bc includes prescription that σ = 1 at the midpoint
"""

struct chebars
    yf::Array{Float64}
    y::Array{Float64}
    D1::Array{Float64,2}
    De::Array{Float64,2}
    Do::Array{Float64,2}
end

mutable struct aring # an abstract, idealized 1D (in r) ring, with periodic BC's and a single forcing moon
    X₀::Float64  # THE key nondimensional parameter.
    ℓ₀::Float64  # The other non-dimensional parameter: sets the size of the domain (in units of Lcrit).
    #
    N::Integer   # number of nodes for a full ring
    nn::Integer  # number of nodes for a half-ring (N÷2 + N%2)
    #
    #ξ::Array{Float64}  # position x normalized to Lcrit
    #Ξ::Array{Float64}  # position x normalized to Lvisc
    #
    β::Float64 # exponent for σ-dependence of viscosity (and diffusivity)
    γ::Float64 # exponent for σ-dependence of relaxation rate 1/s
    #
    chArs::chebars # Tuple of chebdif arrays
    #
    w::Array{Float64}
    σ::Array{Float64}
    Ar::Any
    a::Any
    Br::Any
    b::Any
end

"""
Functions for getting "A" matrix and "a" vector.
"""
function getA(De,σ,w,β,γ)#;hrel=0.0)
    A  = diagm((1+β)*σ.^(β+γ) .- w.^2) * De
    Ar = A[:,2:end-1]
    a  = A[:,end]
    return (Ar,a)
end
function getA!(r::aring)
    (r.Ar,r.a) = getA( r.chArs.De , r.σ , r.w , r.β , r.γ )
end

"""
Functions for getting "B" matrix and "b" vector.
"""
function getB(Do,σ,w,γ,ℓ₀,τ)
    # Get B1.
    B1 = diagm(w.*σ.^γ) ./ (ℓ₀-1) ./ τ
    B1r = B1[:,2:end-1]
    b1  = B1[:,end]
    # And B2.
    B2 = 2 * diagm(w) * diagm(Do*w)
    B2r = B2[:,2:end-1]
    b2  = B2[:,end]
    # Combine
    Br = B1r + B2r
    b  = b1  + b2
    # And return:
    return (Br,b)
end
function getB!(ar::aring)
    τ = ar.X₀^(-1/8)
    (ar.Br, ar.b) = getB( ar.chArs.Do , ar.σ , ar.w , ar.γ , ar.ℓ₀, τ)
end

"""
Functions for setting the normalized velocity "w".
    hrel is the Hill radius relative to Lcrit.
"""
function getw(y,ℓ₀ ; hrel=0.0) # -1 ≤ y ≤ 1
    if ASSERT
        @assert -1 ≤ y ≤ 1
        @assert hrel ≥ 0.0
    end
    ξL = y * (ℓ₀-1) + ℓ₀ #  1 ≤ ξL = position relative to mass on left  (dim'less)
    ξR = y * (ℓ₀-1) - ℓ₀ # -1 ≥ ξR = position relative to mass on right (dim'less)
    # In principle, we should consider a mass at ..., -4L₀, -2L₀, 0, 2L₀, 4L₀, ...,
    #    but honestly, that's overkill in all situations of interest I can think of here.
    # We could also test here whether |ξL| > 2.5*hrel and |ξR| > 2.5*hrel, but we don't, b/c
    #   we might want to violate that in the early stages of an iterative solution.
    # That condition should be tested elswhere, however.
    if hrel ≈ 0.0
        return ξL^-4 - ξR^-4
    else
        aa = 0.711557; bb = -7.58607 # Constants given by Grätz, Seiß & Spahn (2018), p. 3.
        function g(hox)
            min(1/2.5,abs(hox))
        end
        function f(hox)
            (1+ aa*g(hox) + bb*g(hox)^2)^(-1)
        end
        function ww(ξ)
            sign(ξ)/ξ^4 * f(hrel/ξ)
        end
        wL = ww(ξL)
        wR = ww(ξR)
        return wL + wR
    end
end
function getw!(ar::aring ; hrel::Real=0.0)
    ar.w = getw.(ar.chArs.y,ar.ℓ₀;hrel)
    #ar.w = w[1:ar.nn]
end

function iterate!(ar::aring; hrel=0.0)
    getw!(ar; hrel)
    getA!(ar)
    getB!(ar)
    # Solve:
    nn= ar.nn
    print(size(ar.Ar))
    @assert size(ar.Ar) == (nn,nn-2)
    @assert size(ar.a) == (nn,)
    @assert size(ar.Br) == (nn,nn-2)
    @assert size(ar.b) == (nn,)
    σ = (ar.Ar - ar.Br) \ (ar.b - ar.a)
    ar.σ = [0, σ..., 1]
    #ar.σ = σ
    #
    return 0
end

"""
Normal function to initialize an instance of "aring".
"""
function aring(X₀,ℓ₀,N,β,γ; hrel=0.0) # h is the Hill radius in units of Lcrit
    if ASSERT
        @assert X₀ > 2^24           # 2^24 is critical value for X; below this, the method does not apply
        @assert ℓ₀ > 3              # L0/Lcrit needs to be greater than a few; this is a bit of a fuzzy boundary, but certainly > 3.
        @assert ℓ₀ * X₀^(1/24) > 3  # Same for L0/Lvisc; this makes the above line redundant, I know - kept as a mnemonic
    end
    #
    # Find nn, which is half of N, since we treat only half the domain.
    # (Why not just specify nn to start? So we can test/try both odd and even N for the full domain.)
    nn = N÷2 + N%2 # nn: "non-negative" (positions)
    #
    # Get collocation points and differentiation matrices
    # Note that y, D1, De, and Do operate on the domain -1 ≤ y ≤ 1. Need to normalize to -ℓ₀ ≤ ξ ≤ ℓ₀.
    (yf,DM,DeM,DoM) = chebdif(N,1)
    D1 = DM[:,:,1]   # full derivative matrix (N x N)
    y = yf[1:nn]
    De = DeM[:,:,1]   # derivative of even functions (nn x nn)
    Do = DoM[:,:,1]   # derivative of odd functions (nn x nn)
    chArs = chebars(yf,y,D1,De,Do) # for
    #
    # Get (dim'less) velocity:
    w = getw.(y,ℓ₀; hrel)
    #
    # Get "A" and "a"
    σ = ones(nn); # just to start
    (Ar,a) = getA(De,σ,w,β,γ)
    #
    # Get "B" and "b"
    τ = X₀^(-1/8)
    (Br,b) = getB(Do,σ,w,β,τ,ℓ₀)
    #
    ar = aring(X₀, ℓ₀, N, nn, β, γ, chArs, w, σ, Ar, a, Br, b)
    iterate!(ar)
    return ar
end


struct ring
    ar::aring
    moon::Tuple # mass, a0
    disk::Tuple # nu, srel
    X₀::Float64 # fundamental parameter
    β::Float64
    γ::Float64
    Lset::Tuple # lengthscales
    gap::Tuple # width, Δwidth
    h::Float64  # Hill radius
end
function ring(m,β,γ;N=1000,a0=1.0e10,ν=1.0,srel=30,ξ₀=10)#,h=0.0)
    function Ω(a₀)
        return √(G*M♄/a₀^3)
    end
    a₀=a0
    Ω = Ω(a₀)
    α = A₁^2 / 18π * Ω * (m/M♄)^2 * a₀^5
    s₀ = srel / Ω
    κ₀ = 3ν
    X₀ = α^2/κ₀^5/s₀^3
    h = a0 * (m/3/M♄)^(1/3) # Hill radius (in physical units, i.e. cm).

    Lvisc = (α/κ₀)^(1/3)
    w₀ = √(κ₀/s₀)
    Lcrit = (α/w₀)^(1/4)
    Lαs = (α*s₀)^0.2          # another lengthscale constructed from α, ν, s
    Lκs = √(κ₀*s₀)
    Lset = (Lvisc,Lcrit,Lαs,Lκs)

    ar = aring(X₀,ξ₀,N,β,γ; hrel=h/Lcrit)

    for i in 1:5
        iterate!(r)
    end
    # Smooth:
    r.σ = conv(r.σ,[1,4,6,4,1]./16)[3:end-2]
    #

    i10 = findfirst(x->x>0.1,r.σ)
    ξ10 = r.ξ[i10]

    i50 = findfirst(x->x>0.5,r.σ)
    ξ50 = r.ξ[i50]

    i90 = findfirst(x->x>0.9,r.σ)
    ξ90 = r.ξ[i90]

    width = 2*(ξ50 * Lcrit)
    Δwidth = 2*(ξ90-ξ10)*Lcrit
    #width=(0.0,0.0)
    #Δwidth=(0.0,0.0)
    h*=0
    return rring(r,(m,a0),(ν,srel),X₀,β,γ,Lset,(width,Δwidth,2*ξ90*Lcrit,2*ξ10*Lcrit),h)
    # Get 10%, 50%, 90% points
end

end#module
