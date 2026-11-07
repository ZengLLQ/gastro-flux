"""
    Formulation Optimization: Particle Micronization & Dissolution Scaling
"""

using GastroFlux
using Printf

function run_formulation_demo()
    println("="^80)
    println(" GastroFlux: Particle Engineering & Micronization Sensitivity Analysis")
    println("="^80)

    # BCS Class II poorly soluble candidate
    drug = DrugProperties(
        name="Candidate-X",
        molecular_weight=380.0,
        compound_type=:Neutral,
        intrinsic_solubility=0.015, # 15 ug/mL
        peff_jej=2.8e-4,
        density=1.30,
        diffusion_coefficient=5.5e-6
    )

    phys = create_human_physiology(state=:Fasted)
    dose = 100.0 # mg
    radii = [2.0, 5.0, 15.0, 35.0, 75.0] # microns

    println("Simulating oral performance across particle size distributions...")
    println("-"^80)
    @printf("%-16s | %-12s | %-12s | %-12s | %-12s\n",
            "Radius (μm)", "Fa (%)", "Cmax (mg/L)", "Tmax (h)", "AUC (mg*h/L)")
    println("-"^80)

    results = evaluate_particle_size_sensitivity(drug, phys, radii; dose_mg=dose, t_end=24.0)

    for (r, fa, cmax, tmax, auc) in results
        @printf("%14.1f μm | %10.1f%% | %12.3f | %12.2f | %12.2f\n",
                r, fa * 100, cmax, tmax, auc)
    end

    println("="^80)
    println("Biopharmaceutics Insight: Micronization (< 5 μm) elevates dissolution velocity")
    println("via Noyes-Whitney surface area amplification, overcoming solubility-limited absorption.")
    println("="^80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_formulation_demo()
end
