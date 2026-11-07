"""
    ACAT: Advanced Compartmental Absorption & Transit Coupled ODE System
"""

"""
    Index mapping for the 43-dimensional ACAT-PBPK ODE state vector:
    Segments 1:9
    - solid[i]:      index (i - 1) * 4 + 1
    - dissolved[i]:  index (i - 1) * 4 + 2
    - precip[i]:     index (i - 1) * 4 + 3
    - enterocyte[i]: index (i - 1) * 4 + 4

    Systemic PK:
    - 37: central compartment mass (mg)
    - 38: peripheral 1 compartment mass (mg)
    - 39: peripheral 2 compartment mass (mg)
    - 40: cumulative renal eliminated mass (mg)
    - 41: cumulative hepatic eliminated mass (mg)
    - 42: cumulative gut metabolized mass (mg)
    - 43: cumulative portal absorbed mass (mg)
"""
const IDX_SOLID_OFFSET      = 0
const IDX_DISS_OFFSET       = 1
const IDX_PRECIP_OFFSET     = 2
const IDX_ENT_OFFSET        = 3

const IDX_CENTRAL           = 37
const IDX_PERIPH1           = 38
const IDX_PERIPH2           = 39
const IDX_ELIM_RENAL        = 40
const IDX_ELIM_HEPATIC      = 41
const IDX_MET_GUT           = 42
const IDX_ABSORBED_PORTAL   = 43
const IDX_FECAL             = 44

const TOTAL_ODE_VARS        = 44

@inline seg_idx(i::Int, offset::Int) = (i - 1) * 4 + offset + 1

"""
    acat_vector_field(u, p, t)

Evaluates the derivative vector du/dt for the complete ACAT-PBPK model.
`p` is a Tuple: `(drug, form, phys, cl_h, e_h, f_h)`.
"""
function acat_vector_field(u::Vector{Float64}, p, t::Float64)
    drug, form, phys, cl_h, e_h, f_h = p

    du = zeros(Float64, TOTAL_ODE_VARS)
    segs = phys.segments
    n_segs = length(segs)

    # Rate of enterocyte basolateral clearance into mesenteric blood (1/h)
    k_basolateral = 45.0

    total_portal_flux = 0.0

    for i in 1:n_segs
        seg = segs[i]
        m_sol = max(0.0, u[seg_idx(i, IDX_SOLID_OFFSET)])
        m_dis = max(0.0, u[seg_idx(i, IDX_DISS_OFFSET)])
        m_pre = max(0.0, u[seg_idx(i, IDX_PRECIP_OFFSET)])
        m_ent = max(0.0, u[seg_idx(i, IDX_ENT_OFFSET)])

        # Segment volume and approximate enterocyte volume
        vol_lumen = seg.volume_mL
        # Enterocyte volume ~ 5% of luminal volume for small intestine
        vol_ent = max(0.5, vol_lumen * 0.05)

        # 1. Local luminal equilibrium solubility
        s_lum = calculate_solubility(drug, seg.pH, seg.bile_salt_mM, form.bile_micelle_partition)

        # 2. Dissolution rate (mg/h)
        r_diss = dissolution_rate(drug, form, m_sol, m_dis, vol_lumen, s_lum)

        # 3. Precipitation rate (mg/h)
        r_precip = precipitation_rate(form, m_dis, vol_lumen, s_lum)

        # 4. Passive mucosal membrane permeation (mg/h)
        cl_perm = permeability_clearance(drug, seg) # mL/h
        c_lum_dis = m_dis / vol_lumen              # mg/mL
        j_perm = cl_perm * c_lum_dis               # mg/h

        # 5. Active P-gp efflux from enterocyte to lumen (mg/h)
        j_efflux = efflux_rate(drug, seg, m_ent, vol_ent)

        # 6. Enterocyte CYP3A4 metabolic extraction (mg/h)
        r_gut_met = gut_metabolism_rate(drug, seg, m_ent, vol_ent)

        # 7. Enterocyte transfer to portal blood (mg/h)
        j_baso = k_basolateral * m_ent
        total_portal_flux += j_baso

        # Inter-segment transit rate constants
        k_transit = (i == 1) ? (log(2.0) / phys.gastric_half_emptying_h) : (1.0 / seg.transit_time_h)

        # Solid drug balance
        flux_sol_out = k_transit * m_sol
        du[seg_idx(i, IDX_SOLID_OFFSET)] += -r_diss - flux_sol_out

        # Dissolved drug balance
        flux_dis_out = k_transit * m_dis
        du[seg_idx(i, IDX_DISS_OFFSET)] += r_diss - r_precip - j_perm + j_efflux - flux_dis_out

        # Precipitated drug balance
        flux_pre_out = k_transit * m_pre
        du[seg_idx(i, IDX_PRECIP_OFFSET)] += r_precip - flux_pre_out

        # Enterocyte drug balance
        du[seg_idx(i, IDX_ENT_OFFSET)] += j_perm - j_efflux - r_gut_met - j_baso

        # Forward transit into next compartment or unabsorbed fecal excretion
        if i < n_segs
            du[seg_idx(i + 1, IDX_SOLID_OFFSET)] += flux_sol_out
            du[seg_idx(i + 1, IDX_DISS_OFFSET)]  += flux_dis_out
            du[seg_idx(i + 1, IDX_PRECIP_OFFSET)] += flux_pre_out
        else
            du[IDX_FECAL] += flux_sol_out + flux_dis_out + flux_pre_out
        end

        # Track cumulative gut metabolism
        du[IDX_MET_GUT] += r_gut_met
    end

    # Systemic disposition (3-compartment body model)
    m_central = max(0.0, u[IDX_CENTRAL])
    m_p1      = max(0.0, u[IDX_PERIPH1])
    m_p2      = max(0.0, u[IDX_PERIPH2])

    vc = phys.systemic_v_central_L
    c_plasma = m_central / vc # mg/L

    # Portal absorption entering central compartment after hepatic first-pass extraction
    # Fraction Fh escapes liver metabolism: Fh = 1 - Eh
    hepatic_escape_input = f_h * total_portal_flux

    # Systemic clearances
    r_renal = drug.renal_clearance * c_plasma
    r_hepatic_systemic = cl_h * c_plasma

    # Central compartment differential
    du[IDX_CENTRAL] = hepatic_escape_input - (r_renal + r_hepatic_systemic) -
                      phys.systemic_k12_h * m_central + phys.systemic_k21_h * m_p1 -
                      phys.systemic_k13_h * m_central + phys.systemic_k31_h * m_p2

    # Peripheral compartments
    du[IDX_PERIPH1] = phys.systemic_k12_h * m_central - phys.systemic_k21_h * m_p1
    du[IDX_PERIPH2] = phys.systemic_k13_h * m_central - phys.systemic_k31_h * m_p2

    # Cumulative tracker balances
    du[IDX_ELIM_RENAL]       = r_renal
    du[IDX_ELIM_HEPATIC]     = r_hepatic_systemic + (e_h * total_portal_flux)
    du[IDX_ABSORBED_PORTAL]  = total_portal_flux

    return du
end
