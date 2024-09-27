module Interior
"""
I wrote this when, for some reason, I thought I needed to consider thermal
gradients within a chondrule. In fact, the thermal diffusivity of most minerals
is about 1 mm^2/s (e.g. Drury,Geothermics Volume 16, Issue 2, 1987, Pages 105-115),
which means that any mm-sized chondrule will thermally equilibrate in about a
second (!!!) so, everything in here is complete overkill and totally unecessary.

I've kept this here in case I find a use for it at some later time. It's really
just 1D Crank-Nicolson in spherical geometry based on FV formalism.
"""

export T, mote, heat!, cool!, cp, Tμ, speck, coolrad!, evaporate!#, hμ
# models the interior of a chondrule grain (thermal transport, melting, evaporation)
include("Consts.jl")
using .Consts: σSB, ρₛ, cₚ, κ

include("DragAndHeat.jl")
using .DragAndHeat: evap

using LinearAlgebra
const f64 = Float64


J = 1.0e7
kJ = 1.0e3 * J
mol = 24.305 + 28.0855 + 3 * 15.999 # obvs, mol wt may differ a bit bc isotope ratio not same as Earth
ε = 0.9 # emissivity of chondrule (total wild-ass guess)
cp0 = 1.2e7 # will override something?
Δh = 5.0e9 # heat of fusion
Tℓ = 273.15 + 0.5*(1475+1600) # melting pt; just averaging range in Lofgren's old report
h1 = Tℓ * cp0
h2 = Tℓ * cp0 + Δh
ΔH0 = 2812 * kJ / mol # heat of vaporization (estimate)

function cp(T)
    #= specific heat capacity (erg/g/K)

     This is the heat capacity for enstatite (MgSiO₃), from 398K-1000K, from:
        Kenneth M. Krupka
        Bruce S. Hemingway
        Richard A. Robie
        Derrill M. Kerrick
        American Mineralogist (1985) 70 (3-4): 261–271.

    If you look at the curve, it starts to look like nonsense below 150K but it
    at least "looks" reasonable up to 2000K. Not that you can trust it, but you
    can use it for development purposes anyway. In the future, I'll need a better fit.
    =#
    # in J/mol/K
    T = max(min(1000,T),398) # just artificially flattening the curve, for now
    cpJmol = 350.7 - 0.1473*T + 1.679e6*T^(-2) - 4296 * T^(-0.5) + 5.826e-5*T^2
    # cpJmol is about 127.5 J/mol/K at 1000K, which is equivalent to 1.27e7 erg/g/K
    # (note that the molecular weight is 100.3875, so don't be startled by the
    # coincidence in the first 3 significant figures )
    return (cpJmol * J / mol, cpJmol)
end

"""
    Tμ(h): (T,μ)

For working purposes, I need a very crude measure of the temperature T and the
melt fraction μ as a function of specific enthalpy h (or specific internal energy -
same thing, since pressure is basically zero).

I've done a bit of looking into this and the literature is on the one hand dense and
on the other hand hard to generalize from for the specific conditions of interest.

Melt temperature ("liquidus") is 1965C for pure enstatite, but a bit lower for
chondrules. Lofgren, LSPC XXI (old reference), gives 1475-1600C for melt.
Heat of fusion hard to find. Basalt heat of fusion is about 500J/g, equivalent to
about 50kJ/mol (mol of enstatite =100g about). Heating from 0K to 2000K takes,
within a factor of 2x, about 250kJ/mol. So, the heat of fusion isn't huge but it's
not negligible either.

Heat of vaporization IS huge, being about 2812 kJ/mol (minus heat of fusion and heat to get
to whatever temp the thing is). That means that evaporative
cooling will be substantial.

Here, assume a constant cp for solid and liquid of about 1.2e7 erg/g, and assume
heat of fusion of 500J/g = 5.0e9 erg/g.
"""
function Tμ(h)
    if h ≤ h1 # solid
        return (h/cp0,0.0)
    elseif h ≤ h2 # partially melted
        return (Tℓ,(h-h1)/Δh)
    elseif h > h2 # fully melted
        return (Tℓ + (h-h2)/cp0,1.0)
    else # should not get here
        return (NaN,NaN)
    end
end

function hμ(T) # inverse of the above - only needed for initialization
    if T ≤ Tℓ
        return (T*cp0,0.0)
    elseif T > Tℓ
        return (h2 + (T-Tℓ)*cp0,1.0)
    else
        return (NaN,NaN)
    end
end

function T(sie)
    #=
    Not using this yet. Plan is to eventually use a caloric equation of state
    that's a bit more sophisticated than constant cp and E = cp*T.
    =#
    cps = 1.0
    cpℓ = 0.5
    hℓ = 10.0
    he = 30.0
    Tℓ = 1000.0
    Te = 2000.0
    #
    sie0 = cps * Tℓ
    sie1 = sie0 + hℓ
    sie2 = sie1 + cpℓ*(Te-Tℓ)
    sie3 = sie2 + he
    if (sie < sie0)
        return sie / cps
    elseif (sie < sie1)
        return Tℓ
    elseif (sie < sie2)
        return Tℓ + (sie-sie1)/cpℓ
    elseif (sie < sie3)
        return Te
    else
        return NaN # needs some thought...
    end
end

mutable struct speck # physical state of a bulked (unresolved) grain
    D::f64 # Diameter
    m::f64 # mass
    T::f64 # Temperature
    Tmax::f64 # max historical temperature
    E::f64 # (internal) Energy
    μ::f64 # melt fraction
    μmax::f64 # max historical melt fraction
    dead::Bool
end

function speck(D,T)
    m = ρₛ * (π/6) * D^3
    (e,μ) = hμ(T)
    E = e * m
    speck(D,m,T,T,E,μ,μ,false)
end

function heat!(s::speck,ΔE) # Q is a heating *rate*; total energy delivered is Q * Δt
    s.E += ΔE
    (T,μ) = Tμ(s.E/s.m)
    s.T = T
    s.Tmax = max(T,s.Tmax)
    s.μ = μ
    s.μmax = max(μ,s.μmax)
    return nothing
end

heat!(s::speck,Q,Δt) = heat!(s,Q*Δt) # Q is a heating *rate*; total energy delivered is Q * Δt

function coolrad!(s::speck,Δt;Tx=0.0)
    T0 = s.T
    Q = - ε * σSB * π * s.D^2 * (T0^4 - Tx^4)
    heat!(s,Q,Δt)
    T1 = s.T
    # @assert abs((T1-T0)/T0) < 0.02 # if temperature changes by more than 2%, time step too big.
    return nothing
end

#function evaporate!(s::speck,Δm) #Note: a *positive* Δm is a mass *loss*.
#    @assert (Δm < s.m) | (s.m < 1.70e-9) # This is the mass of a 10μm enstatite grain
function evaporate!(s::speck,Δt)
    Δm = evap(s.T) * Δt
    T0 = s.T
    if Δm > s.m
        s.dead = true
        s.m = s.E = s.T = s.D = 0.0
        return nothing
    else
        e1 = s.E / s.m # mean specific energy
        s.m -= Δm
        s.D = (s.m * (6/π) / ρₛ)^(1/3)
        if s.μ ≈ 1.0
            e = e1
        elseif (s.μ > 0.0) & (Δm ≤ s.m * s.μ)
            e = h2
        else
#            @warn "Evaporating solid grain."
            e = e1
        end
        ΔE = Δm * (ΔH0-e) # Heat of vaporiation is relative to STP, so cooling is a bit less effective if you're already hot, obvs
    end
    heat!(s,-ΔE)
    T1 = s.T
#    @assert abs((T1-T0)/T0) < 0.10 # if temperature changes by more than 10%, time step (or mass loss) too big.
end


#function heatadv(s::speck,)
#    Qdot = π * D * κ * ΔT * Nu
#end

mutable struct mote
    D::f64 # Diameter (plan is to update this due to evaporation...)
    m::f64 # mass (plan is to update this due to evaporation...)
    N::Integer # number of cells
    Δt::f64 # suggested time step
    m_::Array{f64,1} # mass interior
    E_::Array{f64,1} # internal energy
    T_::Array{f64,1} # Temperature
    Tmax_::Array{f64,1} # historical max temperature
    Cᵐ::Any #Diagonal{f64,Array{f64,1}} # heat capacity matrix
    Cⁱᵐ::Any# Diagonal{f64,Array{f64,1}} # inverse heat capacity matrix
    Fᵐ::Any#Tridiagonal{f64,Array{f64,1}} # LHS tridiagonal matrix
    #Bᵐ::Any#Tridiagonal{f64,Array{f64,1}} # RHS tridiagonal matrix
    Dᵐ::Any# differencing matrix (is this really needed? probably not...)
    #phase_::Array{f64,1} # phase: 0 = solid, 1 = liquid, 2 = gas
end
function mote(D,N=5,T0=0.0)
    # We could imagine a non-uniform grid, but even though using Crank-Nicholson,
    # time-step is set by von Neumann stability criterion (to ensure accuracy),
    # so there's no real advantage to the added complication.
    # All of the matrices for time-stepping are set here, and "carried" with the
    # struct, to make time-stepping itself easier.
    # NB: R_, M_ = vertex-centered; r_, m_ = cell=centered. (T_, E_ are cell-centered)
    m = ρₛ * (π/6)D^3
    R_ = collect(range(0.0,D/2,length=N+1))[2:end]
    ΔR = (D/2.0)/N
    ΔtVN = 0.5 * ΔR^2 / (κ/ρₛ/cₚ) # von Neumann stability max Δt
    @assert ΔR ≈ R_[2]-R_[1]
    R0_ = [0.0, R_...] # used to construct r_
    r_ = sqrt.((R0_[2:end].^2 + R0_[1:end-1].^2)/2) # too clever by half?
    M_ = m * R_.^3 / (D/2)^3
    m_ = [M_[1],diff(M_)...]
    Cᵐ  = Diagonal(cₚ * m_) # heat capacity matrix
    Cⁱᵐ = Diagonal(1.0./(cₚ * m_)) # inverse heat capacity matrix
    T_ = T0 * ones(N)
    E_ = Cᵐ * T_
    # E_ = e_ .* m_
    # T_ = T.(e_)
    Δt = ΔtVN #* 2
    # constructing k-matrix
    kdiag_ = [R_[1:end-1].^2 ./ diff(r_) ..., 0.0]
    ksuper_ = -kdiag_[1:end-1]
    ksub_ = zeros(N-1) # 0.0 * ksuper_
    kᵐ = 4π * κ * Tridiagonal(ksub_,kdiag_,ksuper_)
    # constructing D-matrix
    Dᵐ = Tridiagonal(ones(N-1),-ones(N),zeros(N-1))
    #
    #α = 0.5
    #Aᵐ = I(N) -    α  * Δt * Dᵐ * kᵐ * Cⁱᵐ
    #Bᵐ = I(N) + (1-α) * Δt * Dᵐ * kᵐ * Cⁱᵐ
    Fᵐ = Dᵐ * kᵐ * Cⁱᵐ
    return mote(D,m,N,Δt,m_,E_,T_,0.0*T_,Cᵐ,Cⁱᵐ,Fᵐ,Dᵐ)
end

"""
    heat!(grain,Q;Δt,α): heats (or cools) grain
    -  Q is the heat power into the grain, integrated over the surface (erg/s)
    -  Δt is the time step (default: 1x von Neumann for grain); can be set arbitrarily
            large - if so, heat!() just iterates until specified Δt is reached. The
            actual Δt in the loop is either 1x the grain's von Neumann Δt, or slightly less.
    - α is weighting for implicit/explicit. α=0 is explicit; α=1 is fully implicit.
"""
function heat!(mot::mote,Q=0.0;Δt=mot.Δt,α=0.5)
    # Δt may be set larger than mote's internal Δt, but not less (i.e., that is overridden).
    # If more, then heat equation is iterated at mote's Δt until specified Δt is reached.
    Ni = ceil(Δt/mot.Δt) # number of iterations
    Δtℓ = Δt / Ni # "loop" Δt
    #E_ = mote.E_;
    Fᵐ = mot.Fᵐ #;  Bᵐ = mot.Bᵐ
    N = mot.N;    Cⁱᵐ = mot.Cⁱᵐ;  Dᵐ = mot.Dᵐ
    #
    Qdot_ = [zeros(N-1)..., -Q]
    P_ = Δtℓ * Dᵐ * Qdot_  # inhomogeneous source vector on RHS
    #= could also use "diff()", in which case:
    Qdot_ = [zeros(N)..., -Q]
    P_ = -Δt * diff(Qdot_) =#
    Aᵐ = I(N) -   α   * Δtℓ * Fᵐ
    Bᵐ = I(N) + (1-α) * Δtℓ * Fᵐ
    E0 = sum(mot.E_)
    for i in 1:Ni
        E_ = mot.E_
        mot.E_ = Aᵐ \ (Bᵐ * E_ + P_)
        mot.T_ = Cⁱᵐ * mot.E_
        mot.Tmax_ = max.(mot.Tmax_,mot.T_)
    end
    E1 = sum(mot.E_)
    ΔE = E1-E0;  ΔQ = Q*Δt
    @assert ΔE ≈ ΔQ
    return nothing #ΔE ≈ ΔQ #nothing
end

function cool!(mot::mote;Δt=mot.Δt,α=0.5)
    Q = -σSB * π * mot.D^2 * mot.T_[end]^4
    heat!(mot,Q;Δt,α)
end

end#module
