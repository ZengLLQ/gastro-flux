using Test
using GastroFlux

@testset "Solubility & Micellar Thermodynamics" begin
    # 1. Weak Acid (e.g. Ibuprofen pKa ~ 4.4, S0 = 0.04 mg/mL)
    acid = DrugProperties(
        name="WeakAcid",
        compound_type=:Acid,
        pKa1=4.4,
        intrinsic_solubility=0.04
    )

    # At pH << pKa, solubility approaches intrinsic solubility S0
    s_low_ph = calculate_solubility(acid, 1.4, 0.0, 0.0)
    @test isapprox(s_low_ph, 0.04, rtol=0.01)

    # At pH == pKa, solubility is 2 * S0 (50% ionized)
    s_pka = calculate_solubility(acid, 4.4, 0.0, 0.0)
    @test isapprox(s_pka, 0.08, rtol=1e-5)

    # At pH = pKa + 2 (pH 6.4), solubility is ~101 * S0
    s_high_ph = calculate_solubility(acid, 6.4, 0.0, 0.0)
    @test isapprox(s_high_ph, 0.04 * (1.0 + 10.0^2.0), rtol=1e-5)

    # 2. Weak Base (e.g. Ketoconazole pKa ~ 6.5, S0 = 0.005 mg/mL)
    base = DrugProperties(
        name="WeakBase",
        compound_type=:Base,
        pKa1=6.5,
        intrinsic_solubility=0.005
    )

    # At pH == pKa, solubility is 2 * S0
    s_base_pka = calculate_solubility(base, 6.5, 0.0, 0.0)
    @test isapprox(s_base_pka, 0.010, rtol=1e-5)

    # In stomach (pH 1.5), base is heavily ionized: S >> S0
    s_base_stomach = calculate_solubility(base, 1.5, 0.0, 0.0)
    @test s_base_stomach > 100.0 * 0.005

    # In small intestine (pH 7.5), base is un-ionized: S approaches S0
    s_base_jej = calculate_solubility(base, 7.5, 0.0, 0.0)
    @test s_base_jej < s_base_stomach

    # 3. Neutral compound
    neutral = DrugProperties(
        name="NeutralDrug",
        compound_type=:Neutral,
        intrinsic_solubility=0.15
    )
    @test calculate_solubility(neutral, 1.5, 0.0, 0.0) == 0.15
    @test calculate_solubility(neutral, 7.4, 0.0, 0.0) == 0.15

    # 4. Bile salt micellar solubilization
    k_mic = 0.08 # mM^-1
    bile_salt_fasted = 3.0 # mM
    bile_salt_fed = 12.0  # mM

    s_aqueous = calculate_solubility(neutral, 6.5, 0.0, k_mic)
    s_mic_fasted = calculate_solubility(neutral, 6.5, bile_salt_fasted, k_mic)
    s_mic_fed = calculate_solubility(neutral, 6.5, bile_salt_fed, k_mic)

    @test s_aqueous == 0.15
    @test s_mic_fasted > s_aqueous
    @test s_mic_fed > s_mic_fasted
    @test isapprox(s_mic_fasted, 0.15 * (1.0 + k_mic * bile_salt_fasted), rtol=1e-5)
    @test isapprox(s_mic_fed, 0.15 * (1.0 + k_mic * bile_salt_fed), rtol=1e-5)
end
