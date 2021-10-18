module rings
using LinearAlgebra, DSP
#include("..\\cheb.jl")
include("cheb.jl")
using .cheb
export ring


"""
Σ is an even function of x
Φ is an odd function of x
U is an odd function of x
"""

# Consts
G  = 6.67430e-8 # Gravitational constant, cgs
M♄ = 5.6834e29  # Saturn, in grams
A₁ = 6.7187     # Goldreich-Tremaine constant
Mm = 1.0e8      # Megameter. Typical orbit is 100Mm.

Pan     = (4.95e19, 1.33584e10) # mass, orbital radius of Pan (Encke gap)
Daphnis = (8.40e16, 1.36504e10) # mass, orbital radius of Daphnis (Keeler gap)

mutable struct ring # an abstract, idealized 1D (in r) ring, with periodic BC's and a single forcing moon
    X::Float64  # THE key nondimensional parameter.
    ξ₀::Float64 # the other key non-dimensional parameter: sets the size of the domain (in dim'less units)
    N::Integer  # number of nodes for a full ring
    β::Float64 # exponent
    nn::Integer # half-ring
    wr::Any
    ξ::Any
    Ξ::Any
    σ::Any
    Ar::Any
    a::Any
    Br::Any
    b::Any
    γ::Float64
end
function ring(X₀,ξ₀,N,β=0.0,γ=0.0;h=0.0)
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
    if h==0.0
        w  = (ξ₀-1)^-4 .* ( (y .+ ξ₀/(ξ₀-1)).^-4 - (y .- ξ₀/(ξ₀-1)).^-4 )
    else
        aa = 0.711557; bb = -7.58607
        xL = y * (ξ₀-1) .+ ξ₀
        xR = y * (ξ₀-1) .- ξ₀
        wL = (xL.^4 + aa*h.*abs.(xL).^3 + bb*h^2*xL.^2).^-1
        wR = (xR.^4 + aa*h.*abs.(xR).^3 + bb*h^2*xR.^2).^-1
        w  = wL - wR
    end
    wr = w[1:nn]
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
    return ring(X₀,ξ₀,N,β,nn,wr,ξξ,ΞΞ,σσ,Ar,a,Br,b,γ)
end

function getA!(r::ring)
    N = r.N; σ = r.σ[2:end]; β = r.β; wr = r.wr; γ=r.γ
    (y,DM,DeM,DoM) = chebdif(N,1)
    D1 = DM[:,:,1]   # full derivative matrix
    De = DeM[:,:,1]  # derivative of even function
    Do = DoM[:,:,1]  # derivative of odd function
    #
    A = diagm((1+β)*σ.^(β+γ) .- wr.^2) * De
    #
    r.Ar = A[:,2:end-1]
    r.a  = A[:,end]
    return 0
end

function getB!(r::ring)
    if r.γ == 0
        return 0
    else
        N = r.N; σ = r.σ[2:end]; β = r.β; wr = r.wr
        (y,DM,DeM,DoM) = chebdif(N,1)
        D1 = DM[:,:,1]   # full derivative matrix
        De = DeM[:,:,1]  # derivative of even function
        Do = DoM[:,:,1]  # derivative of odd function
        #
        τ = r.X₀^(-1/8)
        B1 = diagm(wr.*σ.^γ) ./ (ξ₀-1) ./ τ
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
end
function rring(m,β,γ;N=1000,a0=1.0e10,ν=1.0,srel=30,ξ₀=10,h=0.0)
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
    r = ring(X₀,ξ₀,N,β,γ;h=h)
    Lvisc = (α/κ₀)^(1/3)
    w₀ = √(κ₀/s₀)
    Lcrit = (α/w₀)^(1/4)
    Lαs = (α*s₀)^0.2          # another lengthscale constructed from α, ν, s
    Lκs = √(κ₀*s₀)
    Lset = (Lvisc,Lcrit,Lαs,Lκs)
#    end
#    r = gimmering(m,ν,srel,β,a₀,N,ξ₀=10)
    for i in 1:5
        iterate!(r)
    end
    # Smooth:
    r.σ = conv(r.σ,[1,4,6,4,1]./16)[3:end-2]
    #
    #=
    i10 = findfirst(x->x>0.1,r.σ)
    ξ10 = r.ξ[i10]

    i50 = findfirst(x->x>0.5,r.σ)
    ξ50 = r.ξ[i50]

    i90 = findfirst(x->x>0.9,r.σ)
    ξ90 = r.ξ[i90]

    width = 2*(ξ50 * Lcrit)
    Δwidth = 2*(ξ90-ξ10)*Lcrit
    =#
    width=(0.0,0.0)
    Δwidth=(0.0,0.0)
    return rring(r,(m,a0),(ν,srel),X₀,β,γ,Lset,(width,Δwidth))
    # Get 10%, 50%, 90% points
end

function encke()
    gimmearing(4.95e18,)
end

end#module

#=
function ring(N)
    nn = N÷2 + N%2  # nn: "non-negative" (positions)
    # Get key physical parameters, and lengthscales to non-dimensionalize X-axis.
    α = params[1]; ν = params[2]; s = params[3]; Lrel = params[4] # velocity amplitude, kinematic viscosity, relaxation time, Lrel = L/Lc
    q = 3*ν/s; κ = 3*ν
    NX = α^2/κ^5/s^3
    parmd = Dict("α"=>α, "ν"=>ν, "s"=>s, "q"=>q, "κ"=>κ, "Lrel"=>Lrel, "X"=>NX)
    pa = sort(collect(parmd),by=x->x[2])
    for p in pa
        print("$p\n")
    end
    #
    Lc = (α/√q)^0.25        # position of critical point for hyperbolic transport (could call it "L0" I suppose)
    Lv = (α/κ)^(1.0/3.0) # viscous "critical" point
    Lαs = (α*s)^0.2          # another lengthscale constructed from α, ν, s
    Lνs = √(ν*s)             # another lengthscale      "        "    "   "
    L0 = Lrel*max(Lc,Lv)
    L = L0-Lc
    Ls = Dict("Lcrit"=>Lc, "Lvisc"=>Lv, "Lαs"=>Lαs, "Lνs"=>Lνs, "L"=>L, "L0"=>L0)
    a = sort(collect(Ls),by=x->x[2]) # this is a bit ugly... ugh
    for e in a
        print("$e\n")
    end
    #


    xnn = x[1:nn]
    y = (L+EXTEND*Lc)*x
    ynn = y[1:nn]

    u = (α ./ (y.+L0).^4 ) - (α ./ (y.-L0).^4)
    unn = u[1:nn]

    A = [ Do   0*De  ; 0*Do    De]./(L-Lc)
    AA = A*diagm([unn;unn])
    #DeD = De*diagm()
    B = [ 0*De  -Do  ;  -(κ/s)*De   -(1/s)*diagm(ones(nn))]

    Af = [D1  0*D1 ; 0*D1 D1]./(L-Lc)
    AAf = Af*diagm([u;u])
    Bf =  [0*D1   -D1 ; -(κ/s)*D1   -(1/s)*diagm(ones(N)) ]

    # done.
    return ring(N,nn,parmd,Ls,y,u,A,AA,B,D1,De,Do,Af,AAf,Bf)
end

function setκ!(r::ring,Σ::Array{Float64,1},β::Float64)
    s = r.parmd["s"]; κ = r.parmd["κ"]
    Deκ = r.De * diagm(Σ.^β)
    B = [ 0*r.De  -r.Do  ;  -(κ/s)*Deκ   -(1/s)*diagm(ones(r.nn))]
    r.B = B
end

function Xs(α,κ,s)
    Lc = (α^2 * s / κ)^(1/8)
    Lv = (
end
end#module


=#
