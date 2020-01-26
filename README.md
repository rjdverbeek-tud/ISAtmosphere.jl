# ISAtmosphere

[![Build Status](https://travis-ci.com/rjdverbeek-tud/ISAtmosphere.jl.svg?branch=master)](https://travis-ci.com/rjdverbeek-tud/ISAtmosphere.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rjdverbeek-tud/ISAtmosphere.jl?svg=true)](https://ci.appveyor.com/project/rjdverbeek-tud/ISAtmosphere-jl)
[![Codecov](https://codecov.io/gh/rjdverbeek-tud/ISAtmosphere.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rjdverbeek-tud/ISAtmosphere.jl)
[![Coveralls](https://coveralls.io/repos/github/rjdverbeek-tud/ISAtmosphere.jl/badge.svg?branch=master)](https://coveralls.io/github/rjdverbeek-tud/ISAtmosphere.jl?branch=master)

International Standard Atmospheric (ISA) model

Only metric units are used.

Functions
* T_K     Atmospheric temperature [K]
* p_Pa    Air pressure [Pa]
* ρ_kg_m³ Air density [kg/m³]
* a_m_s   Speed of sound [m/s]
* Vcas2Vtas True airspeed Vtas [m/s] as a function of the calibrated airspeed Vcas [m/s]
* Vtas2Vcas Calibrated airspeed Vcas [m/s] as a function of the True airspeed Vtas [m/s]
* M2Vtas  True airspeed Vtas [m/s] as a function of the Mach number
* Vtas2M  Mach as a function of the True airspeed Vtas [m/s]
* M2Vcas  Calibrated airspeed Vcas [m/s] as a function of the Mach number
* Vtas2M  Mach as a function of the True airspeed Vtas [m/s]
* Hp_trans_m Transition altitude [m]
* conditions return T_K, p_Pa, ρ_kg_m³ and a_m_s in struct AtmosConditions

* θ   temperature ratio
* δ   pressure ratio
* σ   density ratio

Constants
* T₀_K    Standard atmospheric temperature [K] at Mean Sea Level (MSL)
* p₀_Pa   Standard atmospheric pressure [Pa] at Mean Sea Level (MSL)
* ρ₀_kg_m³  Standard atmospheric density [kg/m³] at Mean Sea Level (MSL)
* a₀_m_s  Speed of sound [m/s] at Mean Sea Level (MSL)
* g₀_m_s² Graviation acceleration [m/s²]
* Hp_trop_m Geopotential pressure altitude [m] of Tropopause

Type
* AtmosConditions: T_K, p_Pa, ρ_kg_m³, a_m_s

Source: EUROCONTROL BADA 4 User Manual Chapter 2.2 Atmosphere Model
