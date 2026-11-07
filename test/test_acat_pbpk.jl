using Test
using LinearAlgebra
using GastroFlux

@testset "ACAT & Coupled Systemic PBPK Simulation" begin
    drug = DrugProperties(
        name="ModelAnalgesic",
        molecular_weight=320.0,
        compound_type=:Acid,
        pKa1=4.2,
        intrinsic_solubility=0.08,  # mg/mL
        peff_jej=2.2e-4,            # cm/s
        fraction_unbound_plasma=0.08,
        fraction_unbound_gut=0.9,
        vmax_cyp3a4=2.0,            # mg/h
        km_cyp3a4=0.015,            # mg/mL
        hepatic_clearance_int=60.0, # L/h
        renal_clearance=2.0         # L/h
    )

    dose_val = 150.0 # mg
    form = Formulation(
        dose_mg=dose_val,
        particle_radius_um=20.0,
        precipitation_rate=0.25
    )

    phys = create_human_physiology(state=:Fasted)

    # 1. ODE Vector Field Rigorous Mass Conservation Check
    cl_h, e_h, f_h = hepatic_disposition(drug, phys.liver_blood_flow_L_h)
    p_params = (drug, form, phys, cl_h, e_h, f_h)

    # Test random state vector
    u_test = abs.(randn(GastroFlux.TOTAL_ODE_VARS)) .* 5.0
    du_test = GastroFlux.acat_vector_field(u_test, p_params, 1.0)

    @test length(du_test) == GastroFlux.TOTAL_ODE_VARS

    # The sum of rate of changes of all mass compartments (excluding portal flux index 43)
    # must equal zero analytically by mass balance conservation
    conserved_indices = vcat(1:42, 44)
    total_derivative = sum(du_test[conserved_indices])
    @test isapprox(total_derivative, 0.0, atol=1e-9)

    # 2. Complete Oral Absorption Simulation (RK45)
    res_rk45 = simulate_absorption(drug, form, phys; t_end=24.0, save_dt=0.1)

    @test res_rk45.cmax > 0.0
    @test res_rk45.tmax > 0.0
    @test res_rk45.auc_0_t > 0.0
    @test res_rk45.auc_0_inf >= res_rk45.auc_0_t
    @test 0.0 < res_rk45.fraction_absorbed <= 1.0
    @test 0.0 < res_rk45.bioavailability <= 1.0

    # Test numerical mass conservation across every step of the 24-hour trajectory
    n_steps = length(res_rk45.t)
    for s in 1:n_steps
        mass_sum = sum(res_rk45.u[s, conserved_indices])
        # Integration relative tolerance within 0.05%
        @test isapprox(mass_sum, dose_val, rtol=5e-4)
    end

    # 3. Complete Oral Absorption Simulation with Stiff Rosenbrock Solver
    solver_ros = RosenbrockSolver(reltol=1e-5, abstol=1e-7, max_step=0.1)
    res_ros = simulate_absorption(drug, form, phys; t_end=8.0, save_dt=0.1, solver=solver_ros)

    @test res_ros.cmax > 0.0
    @test res_ros.tmax > 0.0
    @test 0.0 < res_ros.fraction_absorbed <= 1.0
    @test 0.0 < res_ros.bioavailability <= 1.0
    # Bioavailability and Cmax should be qualitatively and quantitatively concordant
    @test isapprox(res_rk45.cmax, res_ros.cmax, rtol=0.20)
    @test isapprox(res_rk45.fraction_absorbed, res_ros.fraction_absorbed, rtol=0.25)
end
