using Test
using Statistics
using GastroFlux

@testset "BCS Classification, NCA & Formulation Sensitivity" begin
    # 1. BCS Class I Model (High Solubility, High Permeability)
    drug_class1 = DrugProperties(
        name="MetoprololLike",
        intrinsic_solubility=5.0, # mg/mL (5000 mg in 250 mL => D0 = 100 / 1250 << 1)
        peff_jej=2.5e-4          # cm/s (High perm: An > 1.0)
    )
    bcs1 = classify_bcs(drug_class1, Formulation(dose_mg=100.0))
    @test bcs1.bcs_class == "Class I"
    @test bcs1.solubility_class == "High"
    @test bcs1.permeability_class == "High"
    @test bcs1.dose_number <= 1.0

    # 2. BCS Class II Model (Low Solubility, High Permeability)
    drug_class2 = DrugProperties(
        name="IbuprofenLike",
        compound_type=:Acid,
        pKa1=4.4,
        intrinsic_solubility=0.005, # mg/mL (S_min at low pH = 0.005 => D0 = 200 / (250*0.005) = 160 >> 1)
        peff_jej=3.0e-4
    )
    bcs2 = classify_bcs(drug_class2, Formulation(dose_mg=200.0))
    @test bcs2.bcs_class == "Class II"
    @test bcs2.solubility_class == "Low"
    @test bcs2.permeability_class == "High"
    @test bcs2.dose_number > 1.0

    # 3. BCS Class III Model (High Solubility, Low Permeability)
    drug_class3 = DrugProperties(
        name="AcyclovirLike",
        intrinsic_solubility=1.5,
        peff_jej=0.2e-4 # Low perm: Peff < 1e-4 cm/s
    )
    bcs3 = classify_bcs(drug_class3, Formulation(dose_mg=100.0))
    @test bcs3.bcs_class == "Class III"
    @test bcs3.solubility_class == "High"
    @test bcs3.permeability_class == "Low"

    # 4. BCS Class IV Model (Low Solubility, Low Permeability)
    drug_class4 = DrugProperties(
        name="FurosemideLike",
        intrinsic_solubility=0.01,
        peff_jej=0.15e-4
    )
    bcs4 = classify_bcs(drug_class4, Formulation(dose_mg=150.0))
    @test bcs4.bcs_class == "Class IV"
    @test bcs4.solubility_class == "Low"
    @test bcs4.permeability_class == "Low"

    # 5. NCA Analytical Benchmark on One-Compartment Disposition: C(t) = C0 * exp(-ke * t)
    c0 = 15.0 # mg/L
    ke = 0.25 # 1/h
    dose = 150.0 # mg
    times = collect(0.0:0.1:24.0)
    c_profile = c0 .* exp.(-ke .* times)

    nca = compute_nca(times, c_profile, dose)

    # Analytical: AUC_inf = C0 / ke = 15 / 0.25 = 60.0
    # MRT = 1 / ke = 4.0 h
    # t_half = ln(2) / ke = 2.77258 h
    @test isapprox(nca.cmax, c0, rtol=1e-4)
    @test isapprox(nca.tmax, 0.0, atol=1e-6)
    @test isapprox(nca.auc_inf, c0 / ke, rtol=1e-3)
    @test isapprox(nca.mrt, 1.0 / ke, rtol=1e-3)
    @test isapprox(nca.t_half, log(2.0) / ke, rtol=1e-3)

    # 6. Particle Size Sensitivity Test
    phys = create_human_physiology(state=:Fasted)
    drug_poor = DrugProperties(
        name="PoorSoluble",
        compound_type=:Neutral,
        intrinsic_solubility=0.02, # Low solubility
        peff_jej=2.5e-4
    )
    sizes = [5.0, 25.0, 100.0] # um
    sens_results = evaluate_particle_size_sensitivity(drug_poor, phys, sizes; dose_mg=100.0, t_end=12.0)

    # Micronized (5 um) should achieve higher Cmax and faster absorption than unmilled (100 um)
    cmax_micronized = sens_results[1][3]
    cmax_unmilled   = sens_results[3][3]
    @test cmax_micronized > cmax_unmilled

    # 7. Food Effect Analysis
    fasted_res, fed_res, cmax_rat, auc_rat = evaluate_food_effect(drug_class1, Formulation(dose_mg=100.0); t_end=12.0)
    # Delayed gastric emptying in fed state leads to later Tmax
    @test fed_res.tmax >= fasted_res.tmax
end
