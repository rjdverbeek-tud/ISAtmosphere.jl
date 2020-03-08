var documenterSearchIndex = {"docs":
[{"location":"#ISAtmosphere.jl-1","page":"Home","title":"ISAtmosphere.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Modules = [ISAtmosphere]\nOrder = [:module]\nPrivate = false","category":"page"},{"location":"#ISAtmosphere.ISAtmosphere","page":"Home","title":"ISAtmosphere.ISAtmosphere","text":"International Standard Atmospheric (ISA) model\n\nThe International Standard Atmosphere (ISA) is a static atmospheric model of how the pressure, temperature, density, and viscosity of the Earth's atmosphere change over a wide range of altitudes or elevations. It has been established to provide a common reference for temperature and pressure and consists of tables of values at various altitudes, plus some formulas by which those values were derived.\n\nIt is also known as the ICAO Standard Atmosphere, ISA is a standard against which to compare the actual atmosphere at any point and time. The real atmosphere differs from ISA in many ways. Sea level pressure varies from day to day, and there are wide extremes of temperature at all levels.\n\nOnly metric units are used.\n\nSource: EUROCONTROL BADA 4 User Manual Chapter 2.2 Atmosphere Model\n\nSource: www.skybrary.aero/index.php/InternationalStandardAtmosphere_(ISA)\n\nSource: en.wikipedia.org/wiki/InternationalStandardAtmosphere\n\n\n\n\n\n","category":"module"},{"location":"#Constants-1","page":"Home","title":"Constants","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Modules = [ISAtmosphere]\nOrder = [:constant]\nPrivate = false","category":"page"},{"location":"#ISAtmosphere.Hp_trop_m","page":"Home","title":"ISAtmosphere.Hp_trop_m","text":"Geopotential pressure altitude [m] of Tropopause\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.R_M²_Ks²","page":"Home","title":"ISAtmosphere.R_M²_Ks²","text":"Real gas constant for air [M²/(Ks²)]\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.T₀_K","page":"Home","title":"ISAtmosphere.T₀_K","text":"Standard atmospheric temperature [K] at Mean Sea Level (MSL)\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.a₀_m_s","page":"Home","title":"ISAtmosphere.a₀_m_s","text":"Speed of sound [m/s] at Mean Sea Level (MSL)\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.g₀_m_s²","page":"Home","title":"ISAtmosphere.g₀_m_s²","text":"Graviation acceleration [m/s²] at Mean Sea Level (MSL)\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.p₀_Pa","page":"Home","title":"ISAtmosphere.p₀_Pa","text":"Standard atmospheric pressure [Pa] at Mean Sea Level (MSL)\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.βT∇_K_m","page":"Home","title":"ISAtmosphere.βT∇_K_m","text":"ISA temperature gradient [K/m] with altitude below the tropopause\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.κ","page":"Home","title":"ISAtmosphere.κ","text":"Adiabatic index of air []\n\n\n\n\n\n","category":"constant"},{"location":"#ISAtmosphere.ρ₀_kg_m³","page":"Home","title":"ISAtmosphere.ρ₀_kg_m³","text":"Standard atmospheric density [kg/m³] at Mean Sea Level (MSL)\n\n\n\n\n\n","category":"constant"},{"location":"#Types-1","page":"Home","title":"Types","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"AtmosConditions","category":"page"},{"location":"#ISAtmosphere.AtmosConditions","page":"Home","title":"ISAtmosphere.AtmosConditions","text":"AtmosConditions(Hp_m, T_K, ΔT_K, p_Pa, ρ_kg_m³, a_m_s)\n\nImmutable STRUCT to keep a set of atmospheric conditions together. This struct can be used to also store an arbitrary set of atmospheric conditions. The function conditions can be used to create the struct.\n\n\n\n\n\n","category":"type"},{"location":"#Functions-1","page":"Home","title":"Functions","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Modules = [ISAtmosphere]\nOrder = [:function]\nPrivate = false","category":"page"},{"location":"#ISAtmosphere.Hp_trans_m-Tuple{Float64,Float64}","page":"Home","title":"ISAtmosphere.Hp_trans_m","text":"Hp_trans_m(Vcas_m_s, M)\n\nReturn the transition altitude Hp_trans_m (also called crossover altitude) [m] between a given calibrated airspeed Vcas_m_s [m/s] and a Mach number M.\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-27/28/29\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.M2Vcas-Tuple{Float64,Float64,Float64}","page":"Home","title":"ISAtmosphere.M2Vcas","text":"M2Vcas(M, p_Pa, T_K)\n\nReturn the calibrated airspeed Vcas_m_s [m/s] as a function of the Mach number M at pressure level p_Pa [Pa] and with temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-26\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.M2Vcas-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.M2Vcas","text":"M2Vcas(M, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.M2Vtas-Tuple{Float64,Float64}","page":"Home","title":"ISAtmosphere.M2Vtas","text":"M2Vtas(M, T_K)\n\nReturn true airspeed Vtas_m_s [m/s] as a function of the Mach number M and temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-26\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.M2Vtas-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.M2Vtas","text":"M2Vtas(M, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.T_K","page":"Home","title":"ISAtmosphere.T_K","text":"T_K(Hp_m [, ΔT_K= 0.0])\n\nReturns the atmospheric temperature [K] at pressure altitude Hp_m [m] and with temperature offset ΔT_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-13/16\n\n\n\n\n\n","category":"function"},{"location":"#ISAtmosphere.Vcas2M-Tuple{Float64,Float64,Float64}","page":"Home","title":"ISAtmosphere.Vcas2M","text":"Vcas2M(Vcas_m_s, T_K)\n\nReturn the Mach number M as a function of the calibrated airspeed Vcas_m_s [m/s] at pressure level p_Pa [Pa] and with temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vcas2M-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.Vcas2M","text":"Vcas2M(Vcas_m_s, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vcas2Vtas-Tuple{Float64,Float64,Float64}","page":"Home","title":"ISAtmosphere.Vcas2Vtas","text":"Vcas2Vtas(Vcas_m_s, p_Pa, T_K)\n\nReturn the true airspeed Vtas_m_s [m/s] as a function of the calibrated airspeed Vcas_m_s [m/s] at pressure level p_Pa [Pa] and with temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-23\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vcas2Vtas-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.Vcas2Vtas","text":"Vcas2Vtas(Vcas_m_s, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vtas2M-Tuple{Float64,Float64}","page":"Home","title":"ISAtmosphere.Vtas2M","text":"Vtas2M(Vtas_m_s, T_K)\n\nReturn the Mach number M as a function of the true airspeed Vtas_m_s [m/s] and temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-26 (reversed)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vtas2M-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.Vtas2M","text":"Vtas2M(Vtas_m_s, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vtas2Vcas-Tuple{Float64,Float64,Float64}","page":"Home","title":"ISAtmosphere.Vtas2Vcas","text":"Vtas2Vcas(Vtas_m_s, p_Pa, T_K)\n\nReturn calibrated airspeed Vcas_m_s [m/s] as a function of the true airspeed Vtas_m_s [m/s] at pressure level p_Pa [Pa] and with temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-24\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.Vtas2Vcas-Tuple{Real,AtmosConditions}","page":"Home","title":"ISAtmosphere.Vtas2Vcas","text":"Vtas2Vcas(Vtas_m_s, AtmosConditions)\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.a_m_s-Tuple{Float64}","page":"Home","title":"ISAtmosphere.a_m_s","text":"a_m_s(T_K)\n\nReturn the speed of sound [m/s] at the temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-22\n\n\n\n\n\n","category":"method"},{"location":"#ISAtmosphere.conditions","page":"Home","title":"ISAtmosphere.conditions","text":"conditions(Hp_m[, ΔT_K = 0.0])\n\nCreate AtmosConditions struct with the atmospheric conditions T_K [K], p_Pa [Pa], and ρ_kg_m³ [kg/m³], and the speed of sound a_m_s [m/s] at a given altitude Hp_m [m] and ΔT_K [K] temperature offset.\n\n\n\n\n\n","category":"function"},{"location":"#ISAtmosphere.p_Pa","page":"Home","title":"ISAtmosphere.p_Pa","text":"p_Pa(Hp_m[, ΔT_K = 0.0])\n\nReturn the air pressure [Pa] at pressure altitude Hp_m [m] and with temperature offset ΔT_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-18/20\n\n\n\n\n\n","category":"function"},{"location":"#ISAtmosphere.ρ_kg_m³-Tuple{Float64,Float64}","page":"Home","title":"ISAtmosphere.ρ_kg_m³","text":"ρ_kg_m³(p_Pa, T_K)\n\nReturn the air density [kg/m³] at pressure level p_Pa [Pa] and temperature T_K [K].\n\nSource: EUROCONTROL BADA 4 User Manual eq. 2.2-21\n\n\n\n\n\n","category":"method"}]
}