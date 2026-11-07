"""
    Sensitivity: Formulation Optimization & Physiological Perturbation Analysis
"""

"""
    evaluate_particle_size_sensitivity(drug::DrugProperties, phys::PhysiologyState, radii_um::Vector{Float64};
                                        dose_mg::Float64 = 100.0, t_end::Float64 = 24.0)

Simulates oral absorption across an array of initial particle sizes (e.g. micronized vs unmilled).
Returns a vector of `(radius_um, fraction_absorbed, cmax, tmax, auc)`.
"""
function evaluate_particle_size_sensitivity(
    drug::DrugProperties,
    phys::PhysiologyState,
    radii_um::Vector{Float64};
    dose_mg::Float64 = 100.0,
    t_end::Float64 = 24.0
)
    results = Tuple{Float64, Float64, Float64, Float64, Float64}[]

    for r in radii_um
        form = Formulation(dose_mg=dose_mg, particle_radius_um=r)
        res = simulate_absorption(drug, form, phys; t_end=t_end)
        push!(results, (r, res.fraction_absorbed, res.cmax, res.tmax, res.auc_0_t))
    end

    return results
end

"""
    evaluate_food_effect(drug::DrugProperties, form::Formulation; t_end::Float64 = 24.0)

Performs parallel simulations under Fasted and Fed physiological conditions.
Returns `(fasted_res::SimulationResult, fed_res::SimulationResult, cmax_ratio::Float64, auc_ratio::Float64)`.
"""
function evaluate_food_effect(
    drug::DrugProperties,
    form::Formulation;
    t_end::Float64 = 24.0
)
    phys_fasted = create_human_physiology(state=:Fasted)
    phys_fed    = create_human_physiology(state=:Fed)

    res_fasted = simulate_absorption(drug, form, phys_fasted; t_end=t_end)
    res_fed    = simulate_absorption(drug, form, phys_fed; t_end=t_end)

    cmax_ratio = res_fed.cmax / max(res_fasted.cmax, 1e-12)
    auc_ratio  = res_fed.auc_0_t / max(res_fasted.auc_0_t, 1e-12)

    return (res_fasted, res_fed, cmax_ratio, auc_ratio)
end
