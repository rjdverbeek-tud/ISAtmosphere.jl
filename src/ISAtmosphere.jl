"""
International Standard Atmospheric (ISA) model

The International Standard Atmosphere (ISA) is a static atmospheric model of how
the pressure, temperature, density, and viscosity of the Earth's atmosphere
change over a wide range of altitudes or elevations. It has been established to
provide a common reference for temperature and pressure and consists of tables
of values at various altitudes, plus some formulas by which those values were
derived.

It is also known as the ICAO Standard Atmosphere, ISA is a standard against
which to compare the actual atmosphere at any point and time. The real
atmosphere differs from ISA in many ways. Sea level pressure varies from day to
day, and there are wide extremes of temperature at all levels.

*Only metric units are used.*

Source: EUROCONTROL BADA 4 User Manual Chapter 2.2 Atmosphere Model

Source: www.skybrary.aero/index.php/International_Standard_Atmosphere_(ISA)

Source: en.wikipedia.org/wiki/International_Standard_Atmosphere
"""
module ISAtmosphere

export T_K, T₀_K, p_Pa, p₀_Pa, ρ_kg_m³, ρ₀_kg_m³, a_m_s, a₀_m_s, g₀_m_s²,
conditions, AtmosConditions, Vcas2Vtas, Vtas2Vcas, M2Vtas, Vtas2M, Hp_trans_m,
M2Vcas, Vcas2M, Hp_trop_m, κ, R_M²_Ks², βT∇_K_m

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
"Graviation acceleration [m/s²] at Mean Sea Level (MSL)"
const g₀_m_s² = 9.80665
"ISA temperature gradient [K/m] with altitude below the tropopause"
const βT∇_K_m = -0.0065
"Geopotential pressure altitude [m] of Tropopause"
const Hp_trop_m = 11000.0

"""
    AtmosConditions(Hp_m, T_K, ΔT_K, p_Pa, ρ_kg_m³, a_m_s)

Immutable STRUCT to keep a set of atmospheric conditions together. This struct
can be used to also store an arbitrary set of atmospheric conditions.
The function `conditions` can be used to create the struct.
"""
struct AtmosConditions
    Hp_m::Float64
    T_K::Float64
    ΔT_K::Float64
    p_Pa::Float64
    ρ_kg_m³::Float64
    a_m_s::Float64
end

"""
    T_K(Hp_m [, ΔT_K= 0.0])

Returns the atmospheric temperature [K] at pressure altitude `Hp_m` [m] and
with temperature offset `ΔT_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-13/16
"""
function T_K(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    if Hp_m ≤ Hp_trop_m
        return T₀_K + ΔT_K + βT∇_K_m * Hp_m
    else
        return T₀_K + ΔT_K + βT∇_K_m * Hp_trop_m
    end
end

T_K(Hp_m::Real, ΔT_K::Real = zero(Real)) = T_K(convert(Float64, Hp_m),
convert(Float64, ΔT_K))

"""
    p_Pa(Hp_m[, ΔT_K = 0.0])

Return the air pressure [Pa] at pressure altitude `Hp_m` [m] and with
temperature offset `ΔT_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-18/20
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

Return the air density [kg/m³] at pressure level `p_Pa` [Pa] and temperature
`T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-21
"""
function ρ_kg_m³(p_Pa::Float64, T_K::Float64)
    return p_Pa / (R_M²_Ks² * T_K)
end

ρ_kg_m³(p_Pa::Real, T_K::Real) = ρ_kg_m³(convert(Float64, p_Pa), convert(Float64, T_K))

"""
    a_m_s(T_K)

Return the speed of sound [m/s] at the temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-22
"""
function a_m_s(T_K::Float64)
    return √(κ*R_M²_Ks²*T_K)
end

a_m_s(T_K::Real) = a_m_s(convert(Float64, T_K))

"""
    Vcas2Vtas(Vcas_m_s, p_Pa, T_K)

Return the true airspeed `Vtas_m_s` [m/s] as a function of the calibrated
airspeed `Vcas_m_s` [m/s] at pressure level `p_Pa` [Pa] and with temperature
`T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-23
"""
function Vcas2Vtas(Vcas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    ρᵥ = ρ_kg_m³(p_Pa, T_K)
    return (2.0*p_Pa/(μ*ρᵥ)*((1.0+p₀_Pa/p_Pa*((1.0+μ*ρ₀_kg_m³/(2.0*p₀_Pa)*Vcas_m_s^2)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vcas2Vtas(Vcas_m_s::Real, p_Pa::Real, T_K::Real) = Vcas2Vtas(convert(Float64, Vcas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
    Vcas2Vtas(Vcas_m_s, AtmosConditions)
"""
Vcas2Vtas(Vcas_m_s::Real, cond::AtmosConditions) = Vcas2Vtas(Vcas_m_s,
cond.p_Pa, cond.T_K)

"""
    Vtas2Vcas(Vtas_m_s, p_Pa, T_K)

Return calibrated airspeed `Vcas_m_s` [m/s] as a function of the true airspeed
`Vtas_m_s` [m/s] at pressure level `p_Pa` [Pa] and with temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-24
"""
function Vtas2Vcas(Vtas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    ρᵥ = ρ_kg_m³(p_Pa, T_K)
    return (2.0*p₀_Pa/(μ*ρ₀_kg_m³)*((1.0+p_Pa/p₀_Pa*((1.0+μ*ρᵥ/(2.0*p_Pa)*Vtas_m_s^2)
            ^(1.0/μ)-1.0))^μ-1.0))^0.5
end

Vtas2Vcas(Vtas_m_s::Real, p_Pa::Real, T_K::Real) = Vtas2Vcas(convert(Float64, Vtas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
    Vtas2Vcas(Vtas_m_s, AtmosConditions)
"""
Vtas2Vcas(Vtas_m_s::Real, conditions::AtmosConditions) = Vtas2Vcas(Vtas_m_s,
conditions.p_Pa, conditions.T_K)

"""
    M2Vtas(M, T_K)

Return true airspeed `Vtas_m_s` [m/s] as a function of the Mach number `M`
and temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vtas(M::Float64, T_K::Float64)
    return M*a_m_s(T_K)
end

M2Vtas(M::Real, T_K::Real) = M2Vtas(convert(Float64, M), convert(Float64, T_K))

M2Vtas(M::Float64, conditions::AtmosConditions) = M * conditions.a_m_s

"""
    M2Vtas(M, AtmosConditions)
"""
M2Vtas(M::Real, conditions::AtmosConditions) = convert(Float64, M) *
conditions.a_m_s

"""
    Vtas2M(Vtas_m_s, T_K)

Return the Mach number `M` as a function of the true airspeed `Vtas_m_s` [m/s]
and temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vtas2M(Vtas_m_s::Float64, T_K::Float64)
    return Vtas_m_s/a_m_s(T_K)
end

Vtas2M(Vtas_m_s::Real, T_K::Real) = Vtas2M(convert(Float64, Vtas_m_s),
convert(Float64, T_K))

"""
    Vtas2M(Vtas_m_s, AtmosConditions)
"""
Vtas2M(Vtas_m_s::Real, conditions::AtmosConditions) = convert(Float64, Vtas_m_s) /
conditions.a_m_s

"""
    M2Vcas(M, p_Pa, T_K)

Return the calibrated airspeed `Vcas_m_s` [m/s] as a function of the Mach number
`M` at pressure level `p_Pa` [Pa] and with temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-26
"""
function M2Vcas(M::Float64, p_Pa::Float64, T_K::Float64)
    return Vtas2Vcas(M2Vtas(M, T_K), p_Pa, T_K)
end

M2Vcas(M::Real, p_Pa::Real, T_K::Real) = M2Vcas(convert(Float64, M),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
    M2Vcas(M, AtmosConditions)
"""
M2Vcas(M::Real, conditions::AtmosConditions) = M2Vcas(M, conditions.p_Pa,
conditions.T_K)

"""
    Vcas2M(Vcas_m_s, T_K)

Return the Mach number `M` as a function of the calibrated airspeed `Vcas_m_s`
[m/s] at pressure level `p_Pa` [Pa] and with temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)
"""
function Vcas2M(Vcas_m_s::Float64, p_Pa::Float64, T_K::Float64)
    return Vtas2M(Vcas2Vtas(Vcas_m_s, p_Pa, T_K), T_K)
end

Vcas2M(Vcas_m_s::Real, p_Pa::Real, T_K::Real) = Vcas2M(convert(Float64, Vcas_m_s),
convert(Float64, p_Pa), convert(Float64, T_K))

"""
    Vcas2M(Vcas_m_s, AtmosConditions)
"""
Vcas2M(Vcas_m_s::Real, conditions::AtmosConditions) = Vcas2M(Vcas_m_s,
conditions.p_Pa, conditions.T_K)

"""
    Hp_trans_m(Vcas_m_s, M)

Return the transition altitude `Hp_trans_m` (also called crossover altitude) [m]
between a given calibrated airspeed `Vcas_m_s` [m/s] and a Mach number `M`.

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-27/28/29
"""
function Hp_trans_m(Vcas_m_s::Float64, M::Float64)
    δ_trans = ((1.0+(κ-1.0)/2.0*(Vcas_m_s/a₀_m_s)^2)
                ^(1.0/μ)-1.0)/((1.0+(κ-1.0)/2.0*M^2)^(1.0/μ)-1.0)
    θ_trans = δ_trans^(-βT∇_K_m * R_M²_Ks² / g₀_m_s²)
    return T₀_K/βT∇_K_m*(θ_trans - 1.0)
end

Hp_trans_m(Vcas_m_s::Real, M::Real) = Hp_trans_m(convert(Float64, Vcas_m_s),
convert(Float64, M))

"""
    θ(T_K)

Return temperature ratio `θ` for temperature `T_K` [K].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-30
"""
θ(T_K::Float64) = T_K / T₀_K
θ(T_K::Real) = θ(convert(Float64, T_K))

"""
    δ(p_Pa)

Return pressure ratio `δ` for air pressure `p_Pa` [Pa].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-31
"""
δ(p_Pa::Float64) = p_Pa / p₀_Pa
δ(p_Pa::Real) = δ(convert(Float64, p_Pa))

"""
    σ(ρ_kg_m³)
Return the density ratio `σ` for the air density `ρ_kg_m³` [kg/m³].

Source: EUROCONTROL BADA 4 User Manual eq. 2.2-32
"""
σ(ρ_kg_m³::Float64) = ρ_kg_m³ / ρ₀_kg_m³
σ(ρ_kg_m³::Real) = σ(convert(Float64, ρ_kg_m³))

"""
    conditions(Hp_m[, ΔT_K = 0.0])

Create `AtmosConditions` struct with the atmospheric conditions `T_K` [K],
`p_Pa` [Pa], and `ρ_kg_m³` [kg/m³], and the speed of sound `a_m_s` [m/s] at a
given altitude `Hp_m` [m] and `ΔT_K` [K] temperature offset.
"""
function conditions(Hp_m::Float64, ΔT_K::Float64 = 0.0)
    T = T_K(Hp_m, ΔT_K)
    p = p_Pa(Hp_m, ΔT_K)
    ρ = ρ_kg_m³(p, T)
    a = a_m_s(T)
    return AtmosConditions(Hp_m, T, ΔT_K, p, ρ, a)
end

end # module