"""
    Dissolution: Noyes-Whitney, Nernst-Brunner & Hixson-Crowell Dynamic Kinetics
"""

"""
    compute_particle_geometry(m_solid::Float64, m_dose::Float64, r0_um::Float64, density::Float64)

Calculates effective particle radius (cm) and total available solid surface area (cm²)
under the Hixson-Crowell cube-root shrinking sphere model:
    r(t) = r_0 * (M(t) / M_0)^(1/3)
    A(t) = 3 * M(t) / (ρ * r(t))
"""
function compute_particle_geometry(
    m_solid::Float64,
    m_dose::Float64,
    r0_um::Float64,
    density::Float64 # g/cm³
)
    if m_solid <= 1e-12 || m_dose <= 1e-12
        return (0.0, 0.0, 0.0)
    end

    r0_cm = r0_um * 1e-4 # um to cm
    density_mg_cm3 = density * 1000.0 # g/cm³ to mg/cm³

    mass_ratio = clamp(m_solid / m_dose, 0.0, 1.0)
    r_current_cm = r0_cm * cbrt(mass_ratio)

    # Avoid division by zero for minuscule radius
    r_effective_cm = max(r_current_cm, 1e-7)

    # Total surface area: A = 3 * M / (rho * r)
    surface_area_cm2 = (3.0 * m_solid) / (density_mg_cm3 * r_effective_cm)

    # Hintz-Johnson boundary layer thickness: δ = min(r, 30 um)
    delta0_cm = 30.0 * 1e-4
    delta_eff_cm = min(r_effective_cm, delta0_cm)

    return (r_effective_cm, surface_area_cm2, delta_eff_cm)
end

"""
    dissolution_rate(drug::DrugProperties, form::Formulation, m_solid::Float64, m_dissolved::Float64, vol_mL::Float64, s_lum::Float64)

Computes instantaneous dissolution mass transfer rate (mg/h) via Nernst-Brunner equation:
    dM/dt = (D / δ) * A * (S_lum - C_diss)
"""
function dissolution_rate(
    drug::DrugProperties,
    form::Formulation,
    m_solid::Float64,
    m_dissolved::Float64,
    vol_mL::Float64,
    s_lum::Float64
)
    if m_solid <= 1e-12 || vol_mL <= 1e-6
        return 0.0
    end

    c_diss = max(0.0, m_dissolved / vol_mL) # mg/mL

    # Driving force: if luminal concentration >= equilibrium solubility, dissolution stops
    driving_force = s_lum - c_diss
    if driving_force <= 0.0
        return 0.0
    end

    r_cm, area_cm2, delta_cm = compute_particle_geometry(
        m_solid, form.dose_mg, form.particle_radius_um, drug.density
    )

    if area_cm2 <= 1e-12 || delta_cm <= 1e-8
        return 0.0
    end

    # D_eff in cm²/h (input is cm²/s)
    d_eff_cm2_h = drug.diffusion_coefficient * 3600.0

    # Rate in mg/h (driving force in mg/mL = mg/cm³)
    rate = (d_eff_cm2_h / delta_cm) * area_cm2 * driving_force

    return max(0.0, rate)
end
