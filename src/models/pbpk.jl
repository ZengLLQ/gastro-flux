"""
    PBPK: Gastrointestinal Absorption & Whole-Body Disposition Orchestrator
"""

"""
    simulate_absorption(drug::DrugProperties, form::Formulation, phys::PhysiologyState;
                        t_end::Float64 = 24.0, save_dt::Float64 = 0.05,
                        solver::AbstractODESolver = RK45Solver(reltol=1e-5, abstol=1e-8))

Runs the complete ACAT coupled biopharmaceutics and 3-compartment systemic disposition simulation.
Returns a `SimulationResult` containing time series and non-compartmental metrics.
"""
function simulate_absorption(
    drug::DrugProperties,
    form::Formulation,
    phys::PhysiologyState;
    t_end::Float64 = 24.0,
    save_dt::Float64 = 0.05,
    solver::AbstractODESolver = RK45Solver(reltol=1e-5, abstol=1e-8)
)
    # Hepatic clearance and extraction
    cl_h, e_h, f_h = hepatic_disposition(drug, phys.liver_blood_flow_L_h)
    p_params = (drug, form, phys, cl_h, e_h, f_h)

    # Initial condition: all dose placed in stomach as solid particles
    u0 = zeros(Float64, TOTAL_ODE_VARS)
    u0[seg_idx(1, IDX_SOLID_OFFSET)] = form.dose_mg

    save_times = collect(0.0:save_dt:t_end)
    if save_times[end] < t_end
        push!(save_times, t_end)
    end

    sol = solve(solver, acat_vector_field, u0, (0.0, t_end); p=p_params, save_times=save_times)

    n_steps = length(sol.t)
    n_segs = length(phys.segments)

    solid_mat = Matrix{Float64}(undef, n_steps, n_segs)
    diss_mat  = Matrix{Float64}(undef, n_steps, n_segs)
    prec_mat  = Matrix{Float64}(undef, n_steps, n_segs)
    ent_mat   = Matrix{Float64}(undef, n_steps, n_segs)

    plasma_conc = Vector{Float64}(undef, n_steps)
    portal_flux = Vector{Float64}(undef, n_steps)

    cum_dissolved = Vector{Float64}(undef, n_steps)
    cum_absorbed  = Vector{Float64}(undef, n_steps)
    cum_gut_met   = Vector{Float64}(undef, n_steps)
    cum_hep_elim  = Vector{Float64}(undef, n_steps)
    cum_ren_elim  = Vector{Float64}(undef, n_steps)
    cum_fecal_elim = Vector{Float64}(undef, n_steps)

    vc = phys.systemic_v_central_L
    k_baso = 45.0

    for step in 1:n_steps
        u_step = sol.u[step]
        flux_sum = 0.0
        diss_sum = 0.0

        for i in 1:n_segs
            sol_val = max(0.0, u_step[seg_idx(i, IDX_SOLID_OFFSET)])
            dis_val = max(0.0, u_step[seg_idx(i, IDX_DISS_OFFSET)])
            pre_val = max(0.0, u_step[seg_idx(i, IDX_PRECIP_OFFSET)])
            ent_val = max(0.0, u_step[seg_idx(i, IDX_ENT_OFFSET)])

            solid_mat[step, i] = sol_val
            diss_mat[step, i]  = dis_val
            prec_mat[step, i]  = pre_val
            ent_mat[step, i]   = ent_val

            diss_sum += dis_val
            flux_sum += k_baso * ent_val
        end

        m_central = max(0.0, u_step[IDX_CENTRAL])
        plasma_conc[step] = m_central / vc # mg/L
        portal_flux[step] = flux_sum

        cum_dissolved[step] = diss_sum
        cum_absorbed[step]  = max(0.0, u_step[IDX_ABSORBED_PORTAL])
        cum_gut_met[step]   = max(0.0, u_step[IDX_MET_GUT])
        cum_hep_elim[step]  = max(0.0, u_step[IDX_ELIM_HEPATIC])
        cum_ren_elim[step]  = max(0.0, u_step[IDX_ELIM_RENAL])
        cum_fecal_elim[step] = max(0.0, u_step[IDX_FECAL])
    end

    # Non-compartmental metrics
    cmax = maximum(plasma_conc)
    max_idx = argmax(plasma_conc)
    tmax = sol.t[max_idx]

    # Trapezoidal AUC 0-t
    auc_0_t = 0.0
    for i in 1:(n_steps - 1)
        dt = sol.t[i + 1] - sol.t[i]
        auc_0_t += 0.5 * (plasma_conc[i] + plasma_conc[i + 1]) * dt
    end

    # Estimate terminal elimination rate constant (kel) from late time points
    start_tail_idx = max(1, floor(Int, 0.75 * n_steps))
    tail_t = sol.t[start_tail_idx:end]
    tail_c = plasma_conc[start_tail_idx:end]

    kel = 0.1
    valid_points = tail_c .> 1e-6
    if count(valid_points) >= 3
        log_c = log.(tail_c[valid_points])
        t_sub = tail_t[valid_points]
        # Linear regression slope: log(C) = intercept - kel * t
        mean_t = mean(t_sub)
        mean_y = mean(log_c)
        denom = sum((t_sub .- mean_t) .^ 2)
        if denom > 1e-12
            slope = sum((t_sub .- mean_t) .* (log_c .- mean_y)) / denom
            if slope < -1e-6
                kel = -slope
            end
        end
    end

    clast = plasma_conc[end]
    auc_extrap = clast / max(kel, 1e-4)
    auc_0_inf = auc_0_t + auc_extrap

    dose = form.dose_mg
    tot_portal = cum_absorbed[end]
    fa = clamp(tot_portal / dose, 0.0, 1.0)

    tot_gut_met = cum_gut_met[end]
    fg = (tot_portal + tot_gut_met > 1e-9) ? clamp(tot_portal / (tot_portal + tot_gut_met), 0.0, 1.0) : 1.0

    f_bio = fa * fg * f_h

    seg_names = [s.name for s in phys.segments]

    return SimulationResult(
        sol.t,
        to_matrix(sol),
        seg_names,
        solid_mat,
        diss_mat,
        prec_mat,
        ent_mat,
        plasma_conc,
        portal_flux,
        cum_dissolved,
        cum_absorbed,
        cum_gut_met,
        cum_hep_elim,
        cum_ren_elim,
        cum_fecal_elim,
        fa,
        fg,
        f_h,
        f_bio,
        cmax,
        tmax,
        auc_0_t,
        auc_0_inf
    )
end
