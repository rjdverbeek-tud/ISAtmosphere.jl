"""
International Standard Atmospheric (ISA) model

Only metric units are used.

Source: EUROCONTROL BADA 4 User Manual Chapter 2.2 Atmosphere Model
"""
module ISAtmosphere

export T, T₀, p, p₀, ρ, ρ₀, a, a₀, g₀, κ, conditions, AtmosConditions,
Vcas2Vtas, Vtas2Vcas, M2Vtas, Vtas2M, Hp_trans, M2Vcas, Vcas2M

"EUROCONTROL BADA 4 User Manual Section 2.2.1/2.2.2"

"Standard atmospheric temperature [K] at Mean Sea Level (MSL)"
const T₀ = 288.15
"Standard atmospheric pressure [Pa] at Mean Sea Level (MSL)"
const p₀ = 101325.0
"Standard atmospheric density [kg/m³] at Mean Sea Level (MSL)"
const ρ₀ = 1.225
"Speed of sound [m/s] at Mean Sea Level (MSL)"
const a₀ = 340.294
"Adiabatic index of air []"
const κ = 1.4
"constant used for Vtas<->Vcas conversion"
const μ = (κ - 1.0)/κ
"Real gas constant for air [M²/(Ks²)]"
const R = 287.05287
"Graviation acceleration [m/s²]"
const g₀ = 9.80665
"ISA temperature gradient [K/m] with altitude below the tropopause"
const βT∇ = -0.0065
"Geopotential pressure altitude [m] of Tropopause"
const Hp_trop = 11000.0

struct AtmosConditions
    T_K::Float64
    p_Pa::Float64
    ρ_kg_m³::Float64
    a_m_s::Float64
end

"""
conditions(Hp_m::Float64, ΔT_K::Float64 = 0.0)

Create struct with the atmospheric conditions T, p, and ρ, and the speed of
sound a at a given altitude [m] and ΔT [K] temperature offset.
"""
function conditions(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    T_K = T(Hp_m, ΔT_K)
    p_Pa = p(Hp_m, ΔT_K)
    ρ_kg_m³ = ρ(p_Pa, T_K)
    a_m_s = a(T_K)
    return AtmosConditions(T_K, p_Pa, ρ_kg_m³, a_m_s)
end

"""
T(Hp::Float64, ΔT::Float64 = 0.0)
Atmospheric temperature [K]
at pressure altitude Hp [m] and with temperature offset ΔT [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-13/16
"""
function T(Hp::Float64, ΔT::Float64 = 0.0)
    if Hp ≤ Hp_trop
        return T₀ + ΔT + βT∇ * Hp
    else
        return T₀ + ΔT + βT∇ * Hp_trop
    end
end

T(Hp::Real, ΔT::Real = zero(Real)) = T(convert(Float64, Hp), convert(Float64, ΔT))

"""
p(Hp::Float64, ΔT::Float64 = 0.0)
Air pressure [Pa]
at pressure altitude Hp [m] and with temperature offset ΔT [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-18/20
"""
function p(Hp::Float64, ΔT::Float64 = 0.0)
    if Hp ≤ Hp_trop
        return p₀*((T(Hp, ΔT) - ΔT) / T₀) ^ (-g₀/(βT∇*R))
    else
        p_trop = p₀*((T(Hp_trop, ΔT) - ΔT) / T₀) ^ (-g₀/(βT∇*R))
        return p_trop * exp(-g₀ / (R*T(Hp_trop)) * (Hp - Hp_trop))
    end
end

p(Hp::Real, ΔT::Real = zero(Real)) = p(convert(Float64, Hp), convert(Float64, ΔT))

"""
ρ(p, T)
Air density [kg/m³] at pressure level [Pa] and Temperature [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-21
"""
function ρ(p::Float64, T::Float64)
    return p / (R * T)
end

ρ(p::Real, T::Real) = ρ(convert(Float64, p), convert(Float64, T))

"""
a(T::Float64)
Speed of sound [m/s]
at pressure level [Pa] and with temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-22
"""
function a(T::Float64)
    return √(κ*R*T)
end

a(T::Real) = a(convert(Float64, T))

"""
Vcas2Vtas(Vcas::Float64, p::Float64, T::Float64)
True airspeed Vtas [m/s] as a function of the calibrated airspeed Vcas [m/s]
at pressure level p [Pa] and with temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-23
"""
function Vcas2Vtas(Vcas::Float64, p::Float64, T::Float64)
    ρᵥ = ρ(p, T)
    return (2.0*p/(μ*ρᵥ)*((1.0+p₀/p*((1.0+μ*ρ₀/(2.0*p₀)*Vcas*Vcas)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vcas2Vtas(Vcas::Real, p::Real, T::Real) = Vcas2Vtas(convert(Float64, Vcas),
convert(Float64, p), convert(Float64, T))

"""
Vtas2Vcas(Vtas::Float64, p::Float64, T::Float64)
Calibrated airspeed Vcas [m/s] as a function of the True airspeed Vtas [m/s]
at pressure level p [Pa] and with temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-24
"""
function Vtas2Vcas(Vtas::Float64, p::Float64, T::Float64)
    ρᵥ = ρ(p, T)
    return (2.0*p₀/(μ*ρ₀)*((1.0+p/p₀*((1.0+μ*ρᵥ/(2.0*p)*Vtas*Vtas)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vtas2Vcas(Vtas::Real, p::Real, T::Real) = Vtas2Vcas(convert(Float64, Vtas),
convert(Float64, p), convert(Float64, T))

"""
M2Vtas(M::Float64, T::Float64)
True airspeed Vtas [m/s] as a function of the Mach number
and temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vtas(M::Float64, T::Float64)
    return M*a(T)
end

M2Vtas(M::Real, T::Real) = M2Vtas(convert(Float64, M), convert(Float64, T))

"""
Vtas2M(Vtas::Float64, T::Float64)
Mach as a function of the True airspeed Vtas [m/s]
and temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vtas2M(Vtas::Float64, T::Float64)
    return Vtas/a(T)
end

Vtas2M(Vtas::Real, T::Real) = Vtas2M(convert(Float64, Vtas), convert(Float64, T))

"""
M2Vcas(M::Float64, p::Float64, T::Float64)
Calibrated airspeed Vcas [m/s] as a function of the Mach number
at pressure level p [Pa] and with temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vcas(M::Float64, p::Float64, T::Float64)
    return Vtas2Vcas(M2Vtas(M, T), p, T)
end

M2Vcas(M::Real, p::Real, T::Real) = M2Vcas(convert(Float64, M),
convert(Float64, p), convert(Float64, T))

"""
Vtas2M(Vtas::Float64, T::Float64)
Mach as a function of the True airspeed Vtas [m/s]
at pressure level p [Pa] and with temperature T [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vcas2M(Vcas::Float64, p::Float64, T::Float64)
    return Vtas2M(Vcas2Vtas(Vcas, p, T), T)
end

Vcas2M(Vcas::Real, p::Real, T::Real) = Vcas2M(convert(Float64, Vcas),
convert(Float64, p), convert(Float64, T))

"""
Hp_trans(Vcas::Float64, M::Float64, ΔT::Float64 = 0.0)
Transition altitude (also called crossover altitude) [m] between
a given calibrated airspeed [m/s] and a Mach number M
and with temperature offset ΔT [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-27/28/29
"""
function Hp_trans(Vcas::Float64, M::Float64, ΔT::Float64 = 0.0)
    δ_trans = ((1.0+(κ-1.0)/2.0*(Vcas/a₀)^2)
                ^(1.0/μ)-1.0)/((1.0+(κ-1.0)/2.0*M^2)^(1.0/μ)-1.0)
    θ_trans = δ_trans^(-βT∇ * R / g₀)
    return T₀/βT∇*(θ_trans - 1.0)
end

Hp_trans(Vcas::Real, M::Real, ΔT::Real = zero(Real)) = Hp_trans(
convert(Float64, Vcas), convert(Float64, M), convert(Float64, ΔT))

"""
θ(T::Float64)
temperature ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-30
"""
θ(T::Float64) = T / T₀
θ(T::Real) = θ(convert(Float64, T))
"""
δ(p::Float64)
pressure ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-31
"""
δ(p::Float64) = p / p₀
δ(p::Real) = δ(convert(Float64, p))
"""
σ(ρ::Float64)
density ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-32
"""
σ(ρ::Float64) = ρ / ρ₀
σ(ρ::Real) = σ(convert(Float64, ρ))

#TODO We recalculate all items all the time. Is there a faster strategy?
#TODO Interpolate Temperature/Pressure data grid

end # module
