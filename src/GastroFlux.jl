module GastroFlux

using LinearAlgebra
using Statistics
using Printf
using Random

# 1. Core Data Structures & Interfaces
include("types.jl")

# 2. High-Accuracy Numerical ODE Solvers
include("solvers/solvers.jl")
include("solvers/rk45.jl")
include("solvers/rosenbrock.jl")

# 3. Biopharmaceutics First-Principles Kinetics
include("biopharmaceutics/solubility.jl")
include("biopharmaceutics/dissolution.jl")
include("biopharmaceutics/precipitation.jl")
include("biopharmaceutics/permeability.jl")

# 4. Gastrointestinal Physiology & Segmental Anatomy
include("physiology/segments.jl")
include("physiology/transporters.jl")
include("physiology/metabolism.jl")

# 5. Coupled ACAT-PBPK Dynamic Models
include("models/acat.jl")
include("models/pbpk.jl")

# 6. Biopharmaceutical Classification & Pharmacokinetic Analysis
include("analysis/bcs.jl")
include("analysis/nca.jl")
include("analysis/sensitivity.jl")

# Exports
export
    # Types
    CompoundType,
    DrugProperties,
    Formulation,
    GISegment,
    PhysiologyState,
    SimulationResult,

    # Solvers
    AbstractODESolver,
    RK45Solver,
    RosenbrockSolver,
    ODESolution,
    solve,
    to_matrix,

    # Biopharmaceutics
    calculate_solubility,
    compute_particle_geometry,
    dissolution_rate,
    precipitation_rate,
    effective_surface_area,
    permeability_clearance,

    # Physiology
    create_human_physiology,
    efflux_rate,
    gut_metabolism_rate,
    hepatic_disposition,

    # Models
    acat_vector_field,
    simulate_absorption,

    # Analysis
    BCSClassification,
    classify_bcs,
    NCAMetrics,
    compute_nca,
    evaluate_particle_size_sensitivity,
    evaluate_food_effect

end # module GastroFlux
