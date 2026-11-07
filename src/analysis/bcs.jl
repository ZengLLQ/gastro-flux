"""
    BCS: Biopharmaceutics Classification System & Dimensionless Transport Numbers
"""

"""
    BCSClassification

Detailed BCS diagnostic result including dimensionless numbers:
- `bcs_class`: "Class I", "Class II", "Class III", or "Class IV"
- `dose_number`: D0 = Dose / (V0 * S_min)
- `dissolution_number`: Dn = (D / δ) * (S_min / (ρ * r0)) * τ_res
- `absorption_number`: An = (Peff / R_lumen) * τ_res
- `solubility_class`: "High" or "Low"
- `permeability_class`: "High" or "Low"
"""
struct BCSClassification
    bcs_class::String
    dose_number::Float64
    dissolution_number::Float64
    absorption_number::Float64
    solubility_class::String
    permeability_class::String
end

"""
    classify_bcs(drug::DrugProperties, form::Formulation; v0_mL::Float64 = 250.0, tau_res_h::Float64 = 3.32)

Computes the microscopic dimensionless numbers (Amidon et al., 1995) and assigns the BCS Class:
- D0 <= 1.0 -> High Solubility
- An >= 1.0 (or Peff >= 1e-4 cm/s) -> High Permeability
"""
function classify_bcs(
    drug::DrugProperties,
    form::Formulation;
    v0_mL::Float64 = 250.0,
    tau_res_h::Float64 = 3.32
)
    # Minimum solubility across GI pH range (1.2 to 7.5)
    test_phs = [1.2, 2.0, 4.5, 6.0, 6.8, 7.5]
    s_min = minimum([calculate_solubility(drug, ph, 0.0, 0.0) for ph in test_phs])

    # 1. Dose Number: D0 = Dose / (V0 * S_min)
    d0 = form.dose_mg / (v0_mL * s_min)

    # 2. Absorption Number: An = (Peff / R_lumen) * tau_res
    # Peff in cm/h, R_lumen ~ 1.0 cm
    peff_cm_h = drug.peff_jej * 3600.0
    r_lumen_cm = 1.0
    an = (peff_cm_h / r_lumen_cm) * tau_res_h

    # 3. Dissolution Number: Dn = (D / delta) * (S_min / (rho * r0)) * tau_res
    d_eff_cm2_h = drug.diffusion_coefficient * 3600.0
    r0_cm = form.particle_radius_um * 1e-4
    delta_cm = min(r0_cm, 30.0 * 1e-4)
    density_mg_cm3 = drug.density * 1000.0

    dn = (d_eff_cm2_h / delta_cm) * (s_min / (density_mg_cm3 * r0_cm)) * tau_res_h

    sol_high = d0 <= 1.0
    perm_high = an >= 1.0

    class_str = if sol_high && perm_high
        "Class I"
    elseif !sol_high && perm_high
        "Class II"
    elseif sol_high && !perm_high
        "Class III"
    else
        "Class IV"
    end

    return BCSClassification(
        class_str,
        d0,
        dn,
        an,
        sol_high ? "High" : "Low",
        perm_high ? "High" : "Low"
    )
end
