"""
    L-Stable 2nd-Order Rosenbrock Stiff ODE Solver (ROS2)
"""

using LinearAlgebra

struct RosenbrockSolver <: AbstractODESolver
    reltol::Float64
    abstol::Float64
    max_steps::Int
    initial_step::Float64
    min_step::Float64
    max_step::Float64
    safety_factor::Float64

    function RosenbrockSolver(;
        reltol::Float64 = 1e-4,
        abstol::Float64 = 1e-7,
        max_steps::Int = 100000,
        initial_step::Float64 = 1e-4,
        min_step::Float64 = 1e-12,
        max_step::Float64 = 0.5,
        safety_factor::Float64 = 0.85
    )
        new(reltol, abstol, max_steps, initial_step, min_step, max_step, safety_factor)
    end
end

"""
    compute_jacobian(f, u, p, t, f0; eps_diff=1e-8)

Evaluates numerical Jacobian J = ∂f/∂u via forward differences.
"""
function compute_jacobian(
    f::Function,
    u::Vector{Float64},
    p,
    t::Float64,
    f0::Vector{Float64};
    eps_diff::Float64 = 1e-7
)
    n = length(u)
    J = Matrix{Float64}(undef, n, n)
    u_pert = copy(u)

    for j in 1:n
        orig_val = u[j]
        h_j = eps_diff * max(abs(orig_val), 1.0)
        u_pert[j] = orig_val + h_j

        f_pert = f(u_pert, p, t)
        for i in 1:n
            J[i, j] = (f_pert[i] - f0[i]) / h_j
        end

        u_pert[j] = orig_val # restore
    end

    return J
end

"""
    solve(solver::RosenbrockSolver, f, u0, tspan; p=nothing, save_times=nothing)

Linearly-implicit L-stable ROS2 solver for stiff differential equations.
"""
function solve(
    solver::RosenbrockSolver,
    f::Function,
    u0::Vector{Float64},
    tspan::Tuple{Float64, Float64};
    p = nothing,
    save_times::Union{Nothing, Vector{Float64}} = nothing
)
    t0, tf = tspan
    @assert tf > t0 "Final time must be greater than initial time"

    n_vars = length(u0)
    t = t0
    u = copy(u0)
    h = min(solver.initial_step, (tf - t0) / 10.0)

    t_history = Float64[t0]
    u_history = Vector{Float64}[copy(u0)]

    # ROS2 Coefficients (L-stable)
    gamma = 1.0 - 1.0 / sqrt(2.0)
    m1 = 1.0 - 1.0 / (2.0 * gamma)
    m2 = 1.0 / (2.0 * gamma)

    n_fevals = 0
    n_steps = 0
    n_rejected = 0

    save_idx = 1
    has_save_times = save_times !== nothing && length(save_times) > 0
    if has_save_times && save_times[1] == t0
        save_idx = 2
    end

    I_mat = Matrix{Float64}(I, n_vars, n_vars)

    while t < tf && n_steps < solver.max_steps
        if t + h > tf
            h = tf - t
        end

        f0 = f(u, p, t)
        n_fevals += 1

        J = compute_jacobian(f, u, p, t, f0)
        n_fevals += n_vars

        # System matrix W = I - gamma * h * J
        W = I_mat - (gamma * h) .* J
        W_lu = lu(W; check=false)

        if !issuccess(W_lu)
            # Jacobian singularity or ill-conditioning, shrink step
            n_rejected += 1
            h = max(h * 0.25, solver.min_step)
            continue
        end

        # Stage 1: W * k1 = f0
        k1 = W_lu \ f0

        # Stage 2: W * k2 = f(u + h*k1, t + h) - 2*gamma*h*J*k1
        u_stage2 = u .+ h .* k1
        f_stage2 = f(u_stage2, p, t + h)
        n_fevals += 1

        rhs2 = f_stage2 .- (2.0 * gamma * h) .* (J * k1)
        k2 = W_lu \ rhs2

        # 2nd-order update
        u_next = u .+ h .* (m1 .* k1 .+ m2 .* k2)

        # Embedded 1st-order error estimate: u_next - (u + h*k1) = h * m2 * (k2 - k1)
        err_vec = h .* m2 .* (k2 .- k1)

        err_norm = 0.0
        for i in 1:n_vars
            sc = solver.abstol + solver.reltol * max(abs(u[i]), abs(u_next[i]))
            err_norm += (err_vec[i] / sc)^2
        end
        err_norm = sqrt(err_norm / n_vars)

        if err_norm <= 1.0 || h <= solver.min_step
            # Step accepted
            t_old = t
            u_old = copy(u)
            t += h
            u = u_next
            n_steps += 1

            if has_save_times
                while save_idx <= length(save_times) && save_times[save_idx] <= t + 1e-12
                    target_t = save_times[save_idx]
                    theta = (target_t - t_old) / (t - t_old)
                    u_interp = (1.0 - theta) .* u_old .+ theta .* u
                    push!(t_history, target_t)
                    push!(u_history, u_interp)
                    save_idx += 1
                end
            else
                push!(t_history, t)
                push!(u_history, copy(u))
            end

            fac = solver.safety_factor * (1.0 / max(err_norm, 1e-4))^0.33
            fac = clamp(fac, 0.2, 3.0)
            h = clamp(h * fac, solver.min_step, solver.max_step)
        else
            # Step rejected
            n_rejected += 1
            fac = solver.safety_factor * (1.0 / err_norm)^0.5
            fac = clamp(fac, 0.1, 0.5)
            h = max(h * fac, solver.min_step)
        end
    end

    retcode = (t >= tf - 1e-10) ? :Success : :MaxStepsExceeded
    return ODESolution(t_history, u_history, retcode, n_steps, n_fevals, n_rejected)
end
