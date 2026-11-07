"""
    Metabolism: Enterocyte CYP3A4 Extraction & Well-Stirred Hepatic Clearance
"""

"""
    gut_metabolism_rate(drug::DrugProperties, seg::GISegment, m_enterocyte::Float64, vol_ent_mL::Float64)

Calculates enterocyte metabolic clearance rate (mg/h) via CYP3A4 Michaelis-Menten kinetics:
    V_met = (Vmax * cyp_rel * C_unbound) / (Km + C_unbound)
"""
function gut_metabolism_rate(
    drug::DrugProperties,
    seg::GISegment,
    m_enterocyte::Float64,
    vol_ent_mL::Float64
)
    if drug.vmax_cyp3a4 <= 0.0 || seg.cyp3a4_relative <= 0.0 || m_enterocyte <= 1e-12 || vol_ent_mL <= 1e-6
        return 0.0
    end

    c_ent = m_enterocyte / vol_ent_mL
    c_unbound = c_ent * drug.fraction_unbound_gut

    vmax_seg = drug.vmax_cyp3a4 * seg.cyp3a4_relative
    rate = (vmax_seg * c_unbound) / (drug.km_cyp3a4 + c_unbound)

    return max(0.0, rate)
end

"""
    hepatic_disposition(drug::DrugProperties, q_liver_L_h::Float64)

Computes hepatic extraction ratio (Eh), hepatic clearance (CL_H, L/h), and
first-pass hepatic escape fraction (Fh) using the well-stirred liver model:
    CL_H = (Q_H * fu,p * CLint / Rb) / (Q_H + fu,p * CLint / Rb)
    E_h = CL_H / Q_H
    F_h = 1 - E_h
"""
function hepatic_disposition(drug::DrugProperties, q_liver_L_h::Float64)
    cl_int = drug.hepatic_clearance_int
    fu_p = drug.fraction_unbound_plasma
    rb = drug.blood_to_plasma_ratio

    eff_cl_int = (fu_p * cl_int) / rb

    cl_h = (q_liver_L_h * eff_cl_int) / (q_liver_L_h + eff_cl_int)
    e_h = cl_h / q_liver_L_h
    f_h = 1.0 - e_h

    return (cl_h, e_h, f_h)
end
