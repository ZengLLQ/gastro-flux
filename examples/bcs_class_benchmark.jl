"""
    BCS Class Benchmark: Comparative Simulation Across Drug Discovery Archetypes
"""

using GastroFlux
using Printf

function run_bcs_benchmark()
    println("="^80)
    println(" GastroFlux: Biopharmaceutical Classification System (BCS) Model Benchmark")
    println("="^80)

    # Class I: High Solubility, High Permeability (e.g., Metoprolol)
    drug_i = DrugProperties(
        name="BCS-I (Metoprolol-like)",
        molecular_weight=267.36,
        compound_type=:Base,
        pKa1=9.7,
        intrinsic_solubility=12.0,   # mg/mL (High)
        peff_jej=2.5e-4,             # cm/s (High)
        hepatic_clearance_int=40.0
    )

    # Class II: Low Solubility, High Permeability (e.g., Ibuprofen)
    drug_ii = DrugProperties(
        name="BCS-II (Ibuprofen-like)",
        molecular_weight=206.29,
        compound_type=:Acid,
        pKa1=4.4,
        intrinsic_solubility=0.04,   # mg/mL (Low)
        peff_jej=3.2e-4,             # cm/s (High)
        hepatic_clearance_int=30.0
    )

    # Class III: High Solubility, Low Permeability (e.g., Acyclovir)
    drug_iii = DrugProperties(
        name="BCS-III (Acyclovir-like)",
        molecular_weight=225.21,
        compound_type=:Ampholyte,
        pKa1=2.27,
        pKa2=9.25,
        intrinsic_solubility=1.6,    # mg/mL (High)
        peff_jej=0.25e-4,            # cm/s (Low)
        hepatic_clearance_int=15.0
    )

    # Class IV: Low Solubility, Low Permeability (e.g., Furosemide)
    drug_iv = DrugProperties(
        name="BCS-IV (Furosemide-like)",
        molecular_weight=330.74,
        compound_type=:Acid,
        pKa1=3.8,
        intrinsic_solubility=0.015,  # mg/mL (Low)
        peff_jej=0.30e-4,            # cm/s (Low)
        hepatic_clearance_int=50.0
    )

    drugs = [drug_i, drug_ii, drug_iii, drug_iv]
    phys = create_human_physiology(state=:Fasted)
    form = Formulation(dose_mg=100.0, particle_radius_um=20.0)

    @printf("%-26s | %-8s | %-6s | %-6s | %-6s | %-6s | %-8s | %-8s\n",
            "Compound", "BCS Class", "Fa (%)", "Fg (%)", "Fh (%)", "F (%)", "Cmax(mg/L)", "AUC(mg*h/L)")
    println("-"^80)

    for d in drugs
        bcs = classify_bcs(d, form)
        res = simulate_absorption(d, form, phys; t_end=24.0)

        @printf("%-26s | %-8s | %6.1f | %6.1f | %6.1f | %6.1f | %10.3f | %10.2f\n",
                d.name, bcs.bcs_class,
                res.fraction_absorbed * 100,
                res.fraction_gut_escape * 100,
                res.fraction_hepatic_escape * 100,
                res.bioavailability * 100,
                res.cmax, res.auc_0_t)
    end
    println("="^80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_bcs_benchmark()
end
