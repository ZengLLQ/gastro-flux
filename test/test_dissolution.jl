using Test
using GastroFlux

@testset "Dissolution & Hixson-Crowell Particle Dynamics" begin
    drug = DrugProperties(
        name="TestCompound",
        intrinsic_solubility=0.1, # mg/mL
        density=1.2,              # g/cm³
        diffusion_coefficient=6.0e-6 # cm²/s
    )

    form = Formulation(
        dose_mg=100.0,
        particle_radius_um=25.0,
        hixson_crowell=true
    )

    # 1. Geometry shrinkage
    r_full, area_full, delta_full = compute_particle_geometry(100.0, 100.0, 25.0, 1.2)
    @test isapprox(r_full, 25.0e-4, rtol=1e-5) # 25 um in cm
    @test area_full > 0.0

    # Half mass remaining
    r_half, area_half, delta_half = compute_particle_geometry(50.0, 100.0, 25.0, 1.2)
    @test isapprox(r_half, 25.0e-4 * cbrt(0.5), rtol=1e-5)
    @test area_half < area_full

    # Zero mass
    r_zero, area_zero, delta_zero = compute_particle_geometry(0.0, 100.0, 25.0, 1.2)
    @test r_zero == 0.0
    @test area_zero == 0.0

    # 2. Dissolution rate under sink conditions
    vol_mL = 250.0
    s_lum = 0.5 # mg/mL
    # Dissolution with 0 dissolved drug
    rate_sink = dissolution_rate(drug, form, 100.0, 0.0, vol_mL, s_lum)
    @test rate_sink > 0.0

    # 3. Saturation: when C_diss >= S_lum, dissolution rate must be zero
    rate_sat = dissolution_rate(drug, form, 100.0, s_lum * vol_mL, vol_mL, s_lum)
    @test rate_sat == 0.0

    # 4. Partial saturation: rate decreases linearly with luminal concentration
    rate_half_sat = dissolution_rate(drug, form, 100.0, 0.5 * s_lum * vol_mL, vol_mL, s_lum)
    @test isapprox(rate_half_sat, 0.5 * rate_sink, rtol=1e-4)

    # 5. Precipitation rate check
    # Undersaturated: no precipitation
    p_rate_under = precipitation_rate(form, 0.5 * s_lum * vol_mL, vol_mL, s_lum)
    @test p_rate_under == 0.0

    # Supersaturated above critical threshold
    m_supersat = 3.0 * s_lum * vol_mL # 3x solubility
    p_rate_super = precipitation_rate(form, m_supersat, vol_mL, s_lum)
    @test p_rate_super > 0.0
end
