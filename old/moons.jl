using Roots

struct moon
    name::String
    a0::Float64 # orbital radius
    Δx::Float64 # gap width
end

m1 = moon("Huygens",117723.5e3,417e3)
m2 = moon("Herschel",118234e3,102e3)
m3 = moon("Russell",118614.5e3,35e3)

moons = (m1,m2,m3)

M♄ = 5.683e26 # kg 
G = 6.6726e-11 # m^3 / s^2 / kg 
A₁ = 6.7187

h(m,a0) = a0 * (m/3/M♄)^(1/3)

Ω(a0) = √(G*M♄/a0^3)

P(a0) = 2π/Ω(a0)

α(m,a0) = A₁^2 / (18π) * Ω(a0) * (m/M♄)^2 * a0^5

function u(x,m,a0)
    xx = abs(x)
    a = 0.71157;  b = -7.58607 
    hh = h(m,a0)
    top = α(m,a0)
    bottom = xx^4 #+ a*hh*xx^3 + b*hh^2*xx^2
    return top/bottom 
end

wavespeed(a0) = √(1.0e-4 / P(a0)) # using low viscosity of 1 cm^2/s

for moon in moons 
    f(m) = u(moon.Δx/2 * 10, m, moon.a0) - wavespeed(moon.a0)
    mm = find_zero(f,(1,1.0e20))
    print(moon.name,"   ",mm)
end
