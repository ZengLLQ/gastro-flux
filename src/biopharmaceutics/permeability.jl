"""
    Permeability: Membrane Mass Transfer & Segmental Surface Area Scaling
"""

"""
    effective_surface_area(seg::GISegment)

Computes the physiologically enlarged mucosal absorptive surface area (cm²)
accounting for anatomical folds (kerckring), villi, and microvillar brush borders:
    A_geom = 2 * π * R * L
    A_abs = A_geom * surface_area_scaling
"""
function effective_surface_area(seg::GISegment)
    if seg.radius_cm <= 0.0 || seg.length_cm <= 0.0
        return 0.0
    end
    geom_area = 2.0 * π * seg.radius_cm * seg.length_cm
    return geom_area * seg.surface_area_scaling
end

"""
    permeability_clearance(drug::DrugProperties, seg::GISegment)

Calculates the passive mucosal absorption clearance (mL/h) in a GI compartment:
    CL_perm = P_eff (cm/h) * A_abs (cm²)
"""
function permeability_clearance(drug::DrugProperties, seg::GISegment)
    peff_cm_s = drug.peff_jej
    peff_cm_h = peff_cm_s * 3600.0

    area_cm2 = effective_surface_area(seg)
    return peff_cm_h * area_cm2
end
