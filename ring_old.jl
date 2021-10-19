module rings
using LinearAlgebra, DSP
#include("..\\cheb.jl")
include("cheb.jl")
using .cheb
export ring


"""

UNITS: cgs

Σ is an even function of x
Φ is an odd function of x
U is an odd function of x
"""

# Consts
G  = 6.67430e-8 # Gravitational constant, cgs
M♄ = 5.6834e29  # mass of Saturn, in grams
A₁ = 6.7187     # Goldreich-Tremaine constant
Mm = 1.0e8      # Megameter. Typical orbit is 100Mm.

Pan     = (4.95e19, 1.33584e10) # (mass, orbital radius) of Pan (Encke gap)
Daphnis = (8.40e16, 1.36504e10) # (mass, orbital radius) of Daphnis (Keeler gap)

mutable struct ring # an abstract, idealized 1D (in r) ring, with periodic BC's and a single forcing moon
    X₀::Float64  # THE key nondimensional parameter.
    ξ₀::Float64 # the other key non-dimensional parameter: sets the size of the domain (in dim'less units)
    N::Integer  # number of nodes for a full ring
    β::Float64 # exponent
    nn::Integer # half-ring
    wr::Any
    wc::Any
    ξ::Any
    Ξ::Any
    σ::Any
    Ar::Any
    a::Any
    Br::Any
    b::Any
    γ::Float64
end
function ring(X₀,ξ₀,N,β=0.0,γ=0.0;h=0.0) # h is the Hill radius in units of Lcrit
    @assert X₀ > 2^24 # 2^24 is critical value for X; below this, the method does not apply
    @assert ξ₀ > 3   # needs to be greater than a few; this is a bit of a fuzzy boundary, but certainly >3, I believe.
    ξvisc = X₀^(1/24)
    #
    # Find nn, which is half of N, since we treat only half the domain.
    # (Why not just specify nn to start? So we can test/try both odd and even N for the full domain.)
    nn = N÷2 + N%2 # nn: "non-negative" (positions)
    #
    # Get collocation points and differentiation matrices
    (y,DM,DeM,DoM) = chebdif(N,1)
    D1 = DM[:,:,1]   # full derivative matrix
    De = DeM[:,:,1]  # derivative of even function
    Do = DoM[:,:,1]  # derivative of odd function
    #
    # Get (dim'less) velocity and reduced (to half domain) velocity
    # Velocity with zero Hill radius:
    w0  = (ξ₀-1)^-4 .* ( (y .+ ξ₀/(ξ₀-1)).^-4 - (y .- ξ₀/(ξ₀-1)).^-4 )
    # Velocity including Hill radius effect:
        aa = 0.711557; bb = -7.58607
        function g(xoh)
            max(2.5,abs(xoh))
        end
        function f(xoh)
            (1+ aa/g(xoh) + bb/g(xoh)/g(xoh))^(-1)
        end
        function ww(x)
            sign(x)/x^4 * f(x/h)
        end
        xL = y * (ξ₀-1) .+ ξ₀
        xR = y * (ξ₀-1) .- ξ₀
        wL = ww.(xL)
        wR = ww.(xR)
        ww  = wL + wR
    wr = w0[1:nn]      # w (right)
    wc = 0.0 * (ww-w0)[1:nn] # w "correction"
    #
    # Form full A matrix. Could use sparse matrix for diagm() but not really needed honestly; we are just constructing
    # the matrix, and memory is not an issue.
    A = diagm((1+β) .-wr.^2) * De
    #
    # Get reduced A and column a.
    Ar = A[:,2:end-1]
    a  = A[:,end]
    @assert size(Ar) == (nn,nn-2)
    @assert size(a) == (nn,)
    #
    # Repeat for B1.
    τ = X₀^(-1/8)
    print("Dimless relaxation time: $τ\n")
    B1 = diagm(wr) ./ (ξ₀-1) ./ τ
    B1r = B1[:,2:end-1]
    b1  = B1[:,end]
    @assert size(B1r) == (nn,nn-2)
    @assert size(b1) == (nn,)
    #
    # And B2.
    B2 = 2 * diagm(wr) * diagm(Do*wr)
    B2r = B2[:,2:end-1]
    b2  = B2[:,end]
    @assert size(B2r) == (nn,nn-2)
    @assert size(b2) == (nn,)
    #
    # Now solve.
    Br = B1r + B2r
    b  = b1  + b2
    σ = (Ar - Br) \ (b - a)
    #
    # And return.
    yr = y[1:nn]
    σ = [0, σ..., 1]
    σσ = [0, σ...] # for plotting purposes
    ξ = ξ₀ .- yr * (ξ₀-1)
    ξξ = [0, ξ...] # for plotting purposes
    ΞΞ = ξξ / ξvisc
    return ring(X₀,ξ₀,N,β,nn,wr,wc,ξξ,ΞΞ,σσ,Ar,a,Br,b,γ)
end

function getA!(r::ring;hrel=0.0)
    N = r.N; σ = r.σ[2:end]; β = r.β; wr = r.wr; wc = r.wc; γ=r.γ
    (y,DM,DeM,DoM) = chebdif(N,1)
    D1 = DM[:,:,1]   # full derivative matrix
    De = DeM[:,:,1]  # derivative of even function
    Do = DoM[:,:,1]  # derivative of odd function
    #
    @assert 0.0 ≤ hrel ≤ 1.0
    w = wr #+ hrel*wc
    A = diagm((1+β)*σ.^(β+γ) .- w.^2) * De
    #
    r.Ar = A[:,2:end-1]
    r.a  = A[:,end]
    return 0
end

function getB!(r::ring;hrel=0.0)
    if r.γ == 0
        return 0
    else
        N = r.N; σ = r.σ[2:end]; β = r.β; wr = r.wr; γ=r.γ; ξ₀=r.ξ₀
        (y,DM,DeM,DoM) = chebdif(N,1)
        D1 = DM[:,:,1]   # full derivative matrix
        De = DeM[:,:,1]  # derivative of even function
        Do = DoM[:,:,1]  # derivative of odd function
        #
        @assert 0.0 ≤ hrel ≤ 1.0
        w = wr #+ hrel*r.wc
        τ = r.X₀^(-1/8)
        B1 = diagm(w.*σ.^γ) ./ (ξ₀-1) ./ τ
        B1r = B1[:,2:end-1]
        b1  = B1[:,end]
        #@assert size(B1r) == (nn,nn-2)
        #@assert size(b1) == (nn,)
        #
        # And B2.
        B2 = 2 * diagm(wr) * diagm(Do*wr)
        B2r = B2[:,2:end-1]
        b2  = B2[:,end]
        #@assert size(B2r) == (nn,nn-2)
        #@assert size(b2) == (nn,)
        #
        # Now solve.
        Br = B1r + B2r
        b  = b1  + b2
        r.Br=Br; r.b=b
    end
    return 0
end

function iterate!(r::ring)
    getA!(r)
    getB!(r)
    Br = r.Br;  Ar = r.Ar;  b = r.b;  a = r.a;
    σ = (Ar - Br) \ (b - a)
    r.σ = [0, 0, σ..., 1]
    return 0
end

#(X₀,ξ₀,N,β=1.0)

struct rring
    r::ring
    moon::Tuple # mass, a0
    disk::Tuple # nu, srel
    X₀::Float64 # fundamental parameter
    β::Float64
    γ::Float64
    Lset::Tuple # lengthscales
    gap::Tuple # width, Δwidth
    h::Float64  # Hill radius
end
function rring(m,β,γ;N=1000,a0=1.0e10,ν=1.0,srel=30,ξ₀=10)#,h=0.0)
#    function gimmearing(m,ν=1.0,srel=30,β=2,a₀=100Mm,N=1000,ξ₀=10;γ=0)
        function Ω(a₀)
            return √(G*M♄/a₀^3)
        end
        a₀=a0
        Ω = Ω(a₀)
        α = A₁^2 / 18π * Ω * (m/M♄)^2 * a₀^5
        s₀ = srel / Ω
        κ₀ = 3ν
        X₀ = α^2/κ₀^5/s₀^3
#        return
    h = a0 * (m/3/M♄)^(1/3)

    Lvisc = (α/κ₀)^(1/3)
    w₀ = √(κ₀/s₀)
    Lcrit = (α/w₀)^(1/4)
    Lαs = (α*s₀)^0.2          # another lengthscale constructed from α, ν, s
    Lκs = √(κ₀*s₀)
    Lset = (Lvisc,Lcrit,Lαs,Lκs)
        r = ring(X₀,ξ₀,N,β,γ;h=h/Lcrit)
#    end
#    r = gimmering(m,ν,srel,β,a₀,N,ξ₀=10)
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
