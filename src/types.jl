"""
    GastroFlux: Advanced Compartmental Absorption & Transit (ACAT) Engine Types
"""

"""
    CompoundType

Classification of ionization behavior for small molecule therapeutics.
Options: `:Acid`, `:Base`, `:Neutral`, `:Ampholyte`.
"""
const CompoundType = Symbol

"""
    DrugProperties

Physicochemical, permeation, metabolic, and disposition parameters of a drug candidate.
All units are standardized:
- `molecular_weight`: g/mol
- `pKa1`, `pKa2`: logarithmic acid dissociation constants (pKa2 = 0.0 if monoprotic)
- `logP`: octanol-water partition coefficient of neutral species
- `intrinsic_solubility`: mg/mL (un-ionized solubility S0)
- `density`: g/cm³ (solid crystal density, typically 1.2 - 1.4 g/cm³)
- `diffusion_coefficient`: cm²/s (effective aqueous diffusivity, ~ 5e-6 cm²/s)
- `peff_jej`: cm/s (human effective jejunal permeability, ~ 1e-4 cm/s)
- `fraction_unbound_plasma`: unbound fraction in systemic circulation (0.0 to 1.0)
- `fraction_unbound_gut`: unbound fraction in enterocytes (0.0 to 1.0)
- `blood_to_plasma_ratio`: blood-to-plasma concentration ratio (Rb)
- `vmax_cyp3a4`: maximum metabolic velocity in enterocytes (mg/h)
- `km_cyp3a4`: Michaelis constant for enterocyte CYP3A4 (mg/mL)
- `vmax_pgp`: maximum active efflux transport rate in enterocytes (mg/h)
- `km_pgp`: Michaelis constant for P-glycoprotein efflux (mg/mL)
- `hepatic_clearance_int`: intrinsic hepatic clearance (L/h)
- `renal_clearance`: renal filtration/secretion clearance (L/h)
"""
struct DrugProperties
    name::String
    molecular_weight::Float64
    compound_type::CompoundType
    pKa1::Float64
    pKa2::Float64
    logP::Float64
    intrinsic_solubility::Float64       # mg/mL
    density::Float64                    # g/cm³
    diffusion_coefficient::Float64      # cm²/s
    peff_jej::Float64                   # cm/s
    fraction_unbound_plasma::Float64
    fraction_unbound_gut::Float64
    blood_to_plasma_ratio::Float64
    vmax_cyp3a4::Float64                # mg/h
    km_cyp3a4::Float64                  # mg/mL
    vmax_pgp::Float64                   # mg/h
    km_pgp::Float64                     # mg/mL
    hepatic_clearance_int::Float64      # L/h
    renal_clearance::Float64            # L/h

    function DrugProperties(;
        name::String = "ModelDrug",
        molecular_weight::Float64 = 350.0,
        compound_type::CompoundType = :Acid,
        pKa1::Float64 = 4.5,
        pKa2::Float64 = 0.0,
        logP::Float64 = 2.5,
        intrinsic_solubility::Float64 = 0.05,
        density::Float64 = 1.25,
        diffusion_coefficient::Float64 = 6.5e-6,
        peff_jej::Float64 = 1.8e-4,
        fraction_unbound_plasma::Float64 = 0.1,
        fraction_unbound_gut::Float64 = 0.8,
        blood_to_plasma_ratio::Float64 = 1.0,
        vmax_cyp3a4::Float64 = 0.0,
        km_cyp3a4::Float64 = 0.01,
        vmax_pgp::Float64 = 0.0,
        km_pgp::Float64 = 0.02,
        hepatic_clearance_int::Float64 = 45.0,
        renal_clearance::Float64 = 1.5
    )
        @assert molecular_weight > 0.0 "Molecular weight must be positive"
        @assert intrinsic_solubility > 0.0 "Intrinsic solubility must be positive"
        @assert density > 0.0 "Density must be positive"
        @assert diffusion_coefficient > 0.0 "Diffusion coefficient must be positive"
        @assert peff_jej >= 0.0 "Permeability must be non-negative"
        @assert 0.0 <= fraction_unbound_plasma <= 1.0 "fu,p must be in [0, 1]"
        @assert 0.0 <= fraction_unbound_gut <= 1.0 "fu,gut must be in [0, 1]"
        @assert blood_to_plasma_ratio > 0.0 "Blood-to-plasma ratio must be positive"

        new(
            name, molecular_weight, compound_type, pKa1, pKa2, logP,
            intrinsic_solubility, density, diffusion_coefficient, peff_jej,
            fraction_unbound_plasma, fraction_unbound_gut, blood_to_plasma_ratio,
            vmax_cyp3a4, km_cyp3a4, vmax_pgp, km_pgp,
            hepatic_clearance_int, renal_clearance
        )
    end
end

"""
    Formulation

Characteristics of the drug product formulation.
- `dose_mg`: administered mass in mg
- `particle_radius_um`: initial mean solid particle radius (microns)
- `precipitation_rate`: first-order or power-law precipitation rate constant (1/h)
- `precipitation_exponent`: exponent β for supersaturation driving force
- `critical_supersaturation`: threshold saturation index (C/Seq - 1) triggering precipitation
- `bile_micelle_partition`: binding/solubilization capacity per mM bile salts (1/mM)
- `hixson_crowell`: whether to model particle surface area shrinkage via cube-root kinetics
"""
struct Formulation
    dose_mg::Float64
    particle_radius_um::Float64
    precipitation_rate::Float64
    precipitation_exponent::Float64
    critical_supersaturation::Float64
    bile_micelle_partition::Float64
    hixson_crowell::Bool

    function Formulation(;
        dose_mg::Float64 = 100.0,
        particle_radius_um::Float64 = 25.0,
        precipitation_rate::Float64 = 0.35,
        precipitation_exponent::Float64 = 1.0,
        critical_supersaturation::Float64 = 1.0,
        bile_micelle_partition::Float64 = 0.065,
        hixson_crowell::Bool = true
    )
        @assert dose_mg > 0.0 "Dose must be positive"
        @assert particle_radius_um > 0.0 "Particle radius must be positive"
        @assert precipitation_rate >= 0.0 "Precipitation rate must be non-negative"
        @assert precipitation_exponent >= 0.0 "Precipitation exponent must be non-negative"

        new(
            dose_mg, particle_radius_um, precipitation_rate,
            precipitation_exponent, critical_supersaturation,
            bile_micelle_partition, hixson_crowell
        )
    end
end

"""
    GISegment

Anatomical, physiological, and biochemical properties of a single GI compartment.
"""
struct GISegment
    name::String
    volume_mL::Float64
    transit_time_h::Float64
    pH::Float64
    bile_salt_mM::Float64
    length_cm::Float64
    radius_cm::Float64
    surface_area_scaling::Float64
    cyp3a4_relative::Float64
    pgp_relative::Float64
end

"""
    PhysiologyState

Complete gastrointestinal tract anatomy and systemic physiological parameters.
Can be configured for `:Fasted` or `:Fed` prandial states.
"""
struct PhysiologyState
    is_fasted::Bool
    gastric_half_emptying_h::Float64
    segments::Vector{GISegment}
    liver_blood_flow_L_h::Float64      # Q_H (~90 L/h in healthy adult)
    portal_blood_flow_L_h::Float64     # Q_portal (~72 L/h)
    systemic_v_central_L::Float64      # V_c (~15 L)
    systemic_v_periph1_L::Float64      # V_p1 (~25 L)
    systemic_v_periph2_L::Float64      # V_p2 (~40 L)
    systemic_k12_h::Float64            # inter-compartmental rate k12 (1/h)
    systemic_k21_h::Float64            # inter-compartmental rate k21 (1/h)
    systemic_k13_h::Float64            # inter-compartmental rate k13 (1/h)
    systemic_k31_h::Float64            # inter-compartmental rate k31 (1/h)
end

"""
    SimulationResult

Output trajectories and summarized pharmacokinetic metrics from GastroFlux.
"""
struct SimulationResult
    t::Vector{Float64}
    u::Matrix{Float64}
    segment_names::Vector{String}
    solid_mass::Matrix{Float64}          # (num_steps, num_segments) mg
    dissolved_mass::Matrix{Float64}      # (num_steps, num_segments) mg
    precipitated_mass::Matrix{Float64}   # (num_steps, num_segments) mg
    enterocyte_mass::Matrix{Float64}     # (num_steps, num_segments) mg
    plasma_concentration::Vector{Float64}# mg/L or ug/mL
    portal_vein_flux::Vector{Float64}    # mg/h
    cumulative_dissolved::Vector{Float64}
    cumulative_absorbed::Vector{Float64}
    cumulative_metabolized_gut::Vector{Float64}
    cumulative_hepatic_eliminated::Vector{Float64}
    cumulative_renal_eliminated::Vector{Float64}
    cumulative_fecal_eliminated::Vector{Float64}
    fraction_absorbed::Float64           # Fa (fraction reaching enterocytes)
    fraction_gut_escape::Float64        # Fg (fraction escaping intestinal metabolism)
    fraction_hepatic_escape::Float64    # Fh (fraction escaping liver first-pass)
    bioavailability::Float64             # F = Fa * Fg * Fh
    cmax::Float64                        # peak plasma concentration (mg/L)
    tmax::Float64                        # time of peak concentration (h)
    auc_0_t::Float64                     # area under curve to t_end (mg*h/L)
    auc_0_inf::Float64                   # extrapolated area under curve (mg*h/L)
end
