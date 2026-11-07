"""
    Precipitation: Supersaturation Driving Force & Recrystallization Dynamics
"""

"""
    precipitation_rate(form::Formulation, m_dissolved::Float64, vol_mL::Float64, s_lum::Float64)

Calculates the mass rate of drug precipitation (mg/h) from supersaturated luminal fluid:
    σ = (C_diss - S_lum) / S_lum
    if σ > σ_crit:
        dM_precip/dt = k_precip * (C_diss - S_lum)^β * V_segment
"""
function precipitation_rate(
    form::Formulation,
    m_dissolved::Float64,
    vol_mL::Float64,
    s_lum::Float64
)
    if m_dissolved <= 1e-12 || vol_mL <= 1e-6 || s_lum <= 1e-12
        return 0.0
    end

    c_diss = m_dissolved / vol_mL # mg/mL

    if c_diss <= s_lum
        return 0.0
    end

    supersat_index = (c_diss - s_lum) / s_lum

    if supersat_index <= form.critical_supersaturation
        return 0.0
    end

    excess_conc = c_diss - (s_lum * (1.0 + form.critical_supersaturation))
    if excess_conc <= 0.0
        return 0.0
    end

    # Rate in mg/h
    rate = form.precipitation_rate * (excess_conc ^ form.precipitation_exponent) * vol_mL

    return max(0.0, rate)
end
