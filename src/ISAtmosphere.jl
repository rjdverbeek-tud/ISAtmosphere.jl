"""
International Standard Atmospheric (ISA) model

Only metric units are used.

Source: EUROCONTROL BADA 4 User Manual Chapter 2.2 Atmosphere Model
"""
module ISAtmosphere

export T_K, T₀_K, p_Pa, p₀_Pa, ρ_kg_m³, ρ₀_kg_m³, a_m_s, a₀_m_s, g₀_m_s²,
conditions, AtmosConditions, Vcas2Vtas, Vtas2Vcas, M2Vtas, Vtas2M, Hp_trans_m,
M2Vcas, Vcas2M, Hp_trop_m

"EUROCONTROL BADA 4 User Manual Section 2.2.1/2.2.2"

"Standard atmospheric temperature [K] at Mean Sea Level (MSL)"
const T₀_K = 288.15
"Standard atmospheric pressure [Pa] at Mean Sea Level (MSL)"
const p₀_Pa = 101325.0
"Standard atmospheric density [kg/m³] at Mean Sea Level (MSL)"
const ρ₀_kg_m³ = 1.225
"Speed of sound [m/s] at Mean Sea Level (MSL)"
const a₀_m_s = 340.294
"Adiabatic index of air []"
const κ = 1.4
"constant used for Vtas<->Vcas conversion"
const μ = (κ - 1.0)/κ
"Real gas constant for air [M²/(Ks²)]"
const R_M²_Ks² = 287.05287
"Graviation acceleration [m/s²]"
const g₀_m_s² = 9.80665
"ISA temperature gradient [K/m] with altitude below the tropopause"
const βT∇_K_m = -0.0065
"Geopotential pressure altitude [m] of Tropopause"
const Hp_trop_m = 11000.0

struct AtmosConditions
    T_K::Float64
    p_Pa::Float64
    ρ_kg_m³::Float64
    a_m_s::Float64
end

"""
T_K(Hp_m::Float64, ΔT_K::Float64 = 0.0)
Atmospheric temperature [K]
at pressure altitude Hp_m [m] and with temperature offset ΔT_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-13/16
"""
function T_K(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    if Hp_m ≤ Hp_trop_m
        return T₀_K + ΔT_K + βT∇_K_m * Hp_m
    else
        return T₀_K + ΔT_K + βT∇_K_m * Hp_trop_m
    end
end

T_K(Hp_m::Real, ΔT_K::Real = zero(Real)) = T_K(convert(Float64, Hp_m), convert(Float64, ΔT_K))

"""
p_Pa(Hp_m::Float64, ΔT_K::Float64 = 0.0)
Air pressure [Pa]
at pressure altitude Hp_m [m] and with temperature offset ΔT_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-18/20
"""
function p_Pa(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    if Hp_m ≤ Hp_trop_m
        return p₀_Pa*((T_K(Hp_m, ΔT_K) - ΔT_K) / T₀_K) ^ (-g₀_m_s²/(βT∇_K_m*R_M²_Ks²))
    else
        p_trop = p₀_Pa*((T_K(Hp_trop_m, ΔT_K) - ΔT_K) / T₀_K) ^ (-g₀_m_s²/(βT∇_K_m*R_M²_Ks²))
        return p_trop * exp(-g₀_m_s² / (R_M²_Ks²*T_K(Hp_trop_m)) * (Hp_m - Hp_trop_m))
    end
end

p_Pa(Hp_m::Real, ΔT_K::Real = zero(Real)) = p_Pa(convert(Float64, Hp_m), convert(Float64, ΔT_K))

"""
ρ_kg_m³(p_Pa, T_K)
Air density [kg/m³] at pressure level [Pa] and Temperature [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-21
"""
function ρ_kg_m³(p_Pa::Float64, T_K::Float64)
    return p_Pa / (R_M²_Ks² * T_K)
end

ρ_kg_m³(p_Pa::Real, T_K::Real) = ρ_kg_m³(convert(Float64, p_Pa), convert(Float64, T_K))

"""
a_m_s(T_K::Float64)
Speed of sound [m/s]
at pressure level [Pa] and with temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-22
"""
function a_m_s(T_K::Float64)
    return √(κ*R_M²_Ks²*T_K)
end

a_m_s(T_K::Real) = a_m_s(convert(Float64, T_K))

"""
Vcas2Vtas(Vcas_m_s::Float64, p_Pa::Float64, T_K::Float64)
True airspeed Vtas_m_s [m/s] as a function of the calibrated airspeed Vcas_m_s [m/s]
at pressure level p_Pa [Pa] and with temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-23
"""
function Vcas2Vtas(Vcas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    ρᵥ = ρ_kg_m³(p_Pa, T_K)
    return (2.0*p_Pa/(μ*ρᵥ)*((1.0+p₀_Pa/p_Pa*((1.0+μ*ρ₀_kg_m³/(2.0*p₀_Pa)*Vcas_m_s^2)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vcas2Vtas(Vcas_m_s::Real, p_Pa::Real, T_K::Real) = Vcas2Vtas(convert(Float64, Vcas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
Vtas2Vcas(Vtas_m_s::Float64, p_Pa::Float64, T_K::Float64)
Calibrated airspeed Vcas_m_s [m/s] as a function of the True airspeed Vtas_m_s [m/s]
at pressure level p_Pa [Pa] and with temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-24
"""
function Vtas2Vcas(Vtas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    ρᵥ = ρ_kg_m³(p_Pa, T_K)
    return (2.0*p₀_Pa/(μ*ρ₀_kg_m³)*((1.0+p_Pa/p₀_Pa*((1.0+μ*ρᵥ/(2.0*p_Pa)*Vtas_m_s^2)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vtas2Vcas(Vtas_m_s::Real, p_Pa::Real, T_K::Real) = Vtas2Vcas(convert(Float64, Vtas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
M2Vtas(M::Float64, T_K::Float64)
True airspeed Vtas_m_s [m/s] as a function of the Mach number
and temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vtas(M::Float64, T_K::Float64)
    return M*a_m_s(T_K)
end

M2Vtas(M::Real, T_K::Real) = M2Vtas(convert(Float64, M), convert(Float64, T_K))

"""
Vtas2M(Vtas_m_s::Float64, T_K::Float64)
Mach as a function of the True airspeed Vtas_m_s [m/s]
and temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vtas2M(Vtas_m_s::Float64, T_K::Float64)
    return Vtas_m_s/a_m_s(T_K)
end

Vtas2M(Vtas_m_s::Real, T_K::Real) = Vtas2M(convert(Float64, Vtas_m_s), convert(Float64, T_K))

"""
M2Vcas(M::Float64, p_Pa::Float64, T_K::Float64)
Calibrated airspeed Vcas_m_s [m/s] as a function of the Mach number
at pressure level p_Pa [Pa] and with temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vcas(M::Float64, p_Pa::Float64, T_K::Float64)
    return Vtas2Vcas(M2Vtas(M, T_K), p_Pa, T_K)
end

M2Vcas(M::Real, p_Pa::Real, T_K::Real) = M2Vcas(convert(Float64, M),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
Vtas2M(Vtas_m_s::Float64, T_K::Float64)
Mach as a function of the True airspeed Vtas_m_s [m/s]
at pressure level p_Pa [Pa] and with temperature T_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vcas2M(Vcas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    return Vtas2M(Vcas2Vtas(Vcas_m_s, p_Pa, T_K), T_K)
end

Vcas2M(Vcas_m_s::Real, p_Pa::Real, T_K::Real) = Vcas2M(convert(Float64, Vcas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
Hp_trans_m(Vcas_m_s::Float64, M::Float64, ΔT_K::Float64 = 0.0)
Transition altitude (also called crossover altitude) [m] between
a given calibrated airspeed [m/s] and a Mach number M
and with temperature offset ΔT_K [K]
EUROCONTROL BADA 4 User Manual eq. 2.2-27/28/29
"""
function Hp_trans_m(Vcas_m_s::Float64, M::Float64, ΔT_K::Float64 = 0.0)
    δ_trans = ((1.0+(κ-1.0)/2.0*(Vcas_m_s/a₀_m_s)^2)
                ^(1.0/μ)-1.0)/((1.0+(κ-1.0)/2.0*M^2)^(1.0/μ)-1.0)
    θ_trans = δ_trans^(-βT∇_K_m * R_M²_Ks² / g₀_m_s²)
    return T₀_K/βT∇_K_m*(θ_trans - 1.0)
end

Hp_trans_m(Vcas_m_s::Real, M::Real, ΔT_K::Real = zero(Real)) = Hp_trans_m(
convert(Float64, Vcas_m_s), convert(Float64, M), convert(Float64, ΔT_K))

"""
θ(T_K::Float64)
temperature ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-30
"""
θ(T_K::Float64) = T_K / T₀_K
θ(T_K::Real) = θ(convert(Float64, T_K))
"""
δ(p_Pa::Float64)
pressure ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-31
"""
δ(p_Pa::Float64) = p_Pa / p₀_Pa
δ(p_Pa::Real) = δ(convert(Float64, p_Pa))
"""
σ(ρ_kg_m³::Float64)
density ratio
EUROCONTROL BADA 4 User Manual eq. 2.2-32
"""
σ(ρ_kg_m³::Float64) = ρ_kg_m³ / ρ₀_kg_m³
σ(ρ_kg_m³::Real) = σ(convert(Float64, ρ_kg_m³))

"""
conditions(Hp_m::Float64, ΔT_K::Float64 = 0.0)

Create struct with the atmospheric conditions T_K, p_Pa, and ρ_kg_m³, and the speed of
sound a at a given altitude [m] and ΔT_K [K] temperature offset.
"""
function conditions(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    T = T_K(Hp_m, ΔT_K)
    p = p_Pa(Hp_m, ΔT_K)
    ρ = ρ_kg_m³(p, T)
    a = a_m_s(T)
    return AtmosConditions(T, p, ρ, a)
end

#TODO We recalculate all items all the time. Is there a faster strategy?
#TODO Interpolate Temperature/Pressure data grid

end # module
