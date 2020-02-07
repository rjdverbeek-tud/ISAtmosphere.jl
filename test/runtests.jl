using ISAtmosphere
using Test

@testset "ISAtmosphere.jl" begin
    @test T_K(9000.0, -1.5) ≈ 228.15 atol = 0.01
    @test T_K(13000.0, 0.0) ≈ 216.65 atol = 0.01
    @test T_K(13000) ≈ 216.65 atol = 0.01
    @test p_Pa(9000.0, -1.5) ≈ 30742.5 atol = 0.1
    @test p_Pa(13000.0, 0.0) ≈ 16510.4 atol = 0.1
    @test p_Pa(13000.0, 0) ≈ 16510.4 atol = 0.1
    @test ρ_kg_m³(30742.5, 228.15) ≈ 0.469414 atol = 0.0001
    @test ρ_kg_m³(16510.4, 216.65) ≈ 0.265483 atol = 0.0001
    @test ρ_kg_m³(16510, 217) ≈ 0.265049 atol = 0.0001
    @test a_m_s(228.15) ≈ 302.8 atol = 0.1
    @test a_m_s(216.65) ≈ 295.07 atol = 0.1
    @test a_m_s(217) ≈ 295.308 atol = 0.01

    @test Vcas2Vtas(128.611, p_Pa(9000.0), T_K(9000.0)) ≈ 201.01 atol = 0.1
    @test Vtas2Vcas(201.01, p_Pa(9000.0), T_K(9000.0)) ≈ 128.611 atol = 0.1
    @test Vtas2Vcas(201, p_Pa(9000.0), T_K(9000.0)) ≈ 128.611 atol = 0.1
    @test Vtas2M(201.01, T_K(9000.0)) ≈ 0.661667 atol = 0.001
    @test Vtas2M(201, T_K(9000.0)) ≈ 0.661667 atol = 0.001
    @test M2Vtas(0.661667, T_K(9000.0)) ≈ 201.01 atol = 0.1
    @test Vcas2M(128.611, p_Pa(9000.0), T_K(9000.0)) ≈ 0.661667 atol = 0.001
    @test M2Vcas(0.661667, p_Pa(9000.0), T_K(9000.0)) ≈ 128.611 atol = 0.1

    @test Hp_trans_m(133.756, 0.8) ≈ 11260.287 atol = 0.1
    @test Hp_trans_m(133.756, 0.6) ≈ 6952.80 atol = 0.1

    @test ISAtmosphere.θ(288) ≈ 0.999479 atol = 0.001
    @test ISAtmosphere.δ(101325) ≈ 1.000 atol = 0.001
    @test ISAtmosphere.σ(1) ≈ 0.816327 atol = 0.001

    cond = ISAtmosphere.conditions(9000.0, -1.5)
    @test cond.Hp_m ≈ 9000.0 atol = 0.1
    @test cond.T_K ≈ 228.15 atol = 0.01
    @test cond.p_Pa ≈ 30742.5 atol = 0.1
    @test cond.ρ_kg_m³ ≈ 0.469414 atol = 0.0001
    @test cond.a_m_s ≈ 302.8 atol = 0.1

    cond2 = conditions(9000.0, 0.0)
    @test cond2.Hp_m ≈ 9000.0 atol = 0.1
    @test cond2.T_K ≈ 229.65 atol = 0.01
    @test cond2.p_Pa ≈ 30742.5 atol = 0.1
    @test cond2.ρ_kg_m³ ≈ 0.466348 atol = 0.0001
    @test cond2.a_m_s ≈ 303.793 atol = 0.1

    @test Vcas2Vtas(128.611, cond2) ≈ 201.01 atol = 0.1
    @test Vtas2Vcas(201.01, cond2) ≈ 128.611 atol = 0.1
    @test Vtas2Vcas(201, cond2) ≈ 128.611 atol = 0.1
    @test Vtas2M(201.01, cond2) ≈ 0.661667 atol = 0.001
    @test Vtas2M(201, cond2) ≈ 0.661667 atol = 0.001
    @test M2Vtas(0.661667, cond2) ≈ 201.01 atol = 0.1
    @test Vcas2M(128.611, cond2) ≈ 0.661667 atol = 0.001
    @test M2Vcas(0.661667, cond2) ≈ 128.611 atol = 0.1
end
