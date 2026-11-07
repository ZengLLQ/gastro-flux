"""
    Fed vs. Fasted: Food Effect & Bile Solubilization Dynamics
"""

using GastroFlux
using Printf

function run_food_effect_demo()
    println("="^80)
    println(" GastroFlux: Food Effect & Prandial State Dynamics on Weak Base Bioavailability")
    println("="^80)

    # Lipophilic poorly soluble weak base (exhibiting positive food effect from bile release)
    lipophilic_base = DrugProperties(
        name="Ketoconazole-like",
        molecular_weight=531.43,
        compound_type=:Base,
        pKa1=6.5,
        intrinsic_solubility=0.003, # mg/mL (Extremely lipophilic & sparingly soluble)
        peff_jej=3.5e-4,
        fraction_unbound_plasma=0.05,
        hepatic_clearance_int=70.0
    )

    form = Formulation(
        dose_mg=200.0,
        particle_radius_um=15.0,
        bile_micelle_partition=0.10 # High bile micellar partition coefficient
    )

    fasted_res, fed_res, cmax_ratio, auc_ratio = evaluate_food_effect(lipophilic_base, form; t_end=24.0)

    println("Simulation Results:")
    @printf("  Fasted State: Cmax = %6.3f mg/L, Tmax = %4.2f h, AUC = %6.2f mg*h/L, Fa = %5.1f%%\n",
            fasted_res.cmax, fasted_res.tmax, fasted_res.auc_0_t, fasted_res.fraction_absorbed * 100)
    @printf("  Fed State:    Cmax = %6.3f mg/L, Tmax = %4.2f h, AUC = %6.2f mg*h/L, Fa = %5.1f%%\n",
            fed_res.cmax, fed_res.tmax, fed_res.auc_0_t, fed_res.fraction_absorbed * 100)
    println("-"^80)
    @printf("  Prandial Impact: Cmax Fed/Fasted Ratio = %4.2fx, AUC Fed/Fasted Ratio = %4.2fx\n",
            cmax_ratio, auc_ratio)
    println("  Conclusion: High fat meal triggers gallbladder bile secretion (12 mM vs 3.5 mM),")
    println("  substantially enhancing micellar solubilization and driving a positive food effect.")
    println("="^80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_food_effect_demo()
end
