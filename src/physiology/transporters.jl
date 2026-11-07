"""
    Transporters: Active Carrier Kinetics & Apical Efflux Dynamics (P-gp)
"""

"""
    efflux_rate(drug::DrugProperties, seg::GISegment, m_enterocyte::Float64, vol_ent_mL::Float64)

Calculates active P-glycoprotein (ABCB1) mediated apical efflux rate (mg/h)
from enterocyte cytosol back into the intestinal lumen:
    J_efflux = (Vmax * pgp_rel * C_unbound) / (Km + C_unbound)
"""
function efflux_rate(
    drug::DrugProperties,
    seg::GISegment,
    m_enterocyte::Float64,
    vol_ent_mL::Float64
)
    if drug.vmax_pgp <= 0.0 || seg.pgp_relative <= 0.0 || m_enterocyte <= 1e-12 || vol_ent_mL <= 1e-6
        return 0.0
    end

    c_ent = m_enterocyte / vol_ent_mL # mg/mL
    c_unbound = c_ent * drug.fraction_unbound_gut

    vmax_seg = drug.vmax_pgp * seg.pgp_relative
    rate = (vmax_seg * c_unbound) / (drug.km_pgp + c_unbound)

    return max(0.0, rate)
end
