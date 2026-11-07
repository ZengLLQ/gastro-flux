"""
    Solubility: Ionization Equilibria & Micellar Solubilization Mechanics
"""

"""
    calculate_solubility(drug::DrugProperties, pH::Float64, bile_salt_mM::Float64, k_mic::Float64)

Computes total equilibrium solubility (mg/mL) in a GI compartment at given pH and bile salt level.
Incorporates Henderson-Hasselbalch ionization and micellar partitioning:
    S_aq(pH) = S_0 * α_ion(pH)
    S_lum = S_aq(pH) * (1 + K_mic * [BileSalt])
"""
function calculate_solubility(
    drug::DrugProperties,
    pH::Float64,
    bile_salt_mM::Float64,
    k_mic::Float64
)
    s0 = drug.intrinsic_solubility
    c_type = drug.compound_type

    alpha_ion = 1.0

    if c_type == :Acid
        # Weak acid: S = S0 * (1 + 10^(pH - pKa))
        delta = pH - drug.pKa1
        # Clamp exponent to prevent numerical overflow
        delta_clamped = clamp(delta, -10.0, 10.0)
        alpha_ion = 1.0 + 10.0^delta_clamped

    elseif c_type == :Base
        # Weak base: S = S0 * (1 + 10^(pKa - pH))
        delta = drug.pKa1 - pH
        delta_clamped = clamp(delta, -10.0, 10.0)
        alpha_ion = 1.0 + 10.0^delta_clamped

    elseif c_type == :Ampholyte
        # Zwitterion / Ampholyte with two pKa values (pKa1 < pKa2)
        # S = S0 * (1 + 10^(pKa1 - pH) + 10^(pH - pKa2))
        d1 = clamp(drug.pKa1 - pH, -10.0, 10.0)
        d2 = clamp(pH - drug.pKa2, -10.0, 10.0)
        alpha_ion = 1.0 + 10.0^d1 + 10.0^d2

    elseif c_type == :Neutral
        alpha_ion = 1.0
    end

    s_aqueous = s0 * alpha_ion

    # Micellar enhancement from bile salts
    bile_effect = 1.0 + max(0.0, k_mic * bile_salt_mM)

    return s_aqueous * bile_effect
end
