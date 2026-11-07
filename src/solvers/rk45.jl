"""
    Dormand-Prince 5(4) Adaptive Runge-Kutta ODE Solver (RK45)
"""

struct RK45Solver <: AbstractODESolver
    reltol::Float64
    abstol::Float64
    max_steps::Int
    initial_step::Float64
    min_step::Float64
    max_step::Float64
    safety_factor::Float64

    function RK45Solver(;
        reltol::Float64 = 1e-6,
        abstol::Float64 = 1e-9,
        max_steps::Int = 100000,
        initial_step::Float64 = 1e-3,
        min_step::Float64 = 1e-12,
        max_step::Float64 = 1.0,
        safety_factor::Float64 = 0.9
    )
        new(reltol, abstol, max_steps, initial_step, min_step, max_step, safety_factor)
    end
end

"""
    solve(solver::RK45Solver, f, u0, tspan; save_times=nothing)

Integrates `du/dt = f(u, p, t)` from `tspan[1]` to `tspan[2]` starting at `u0`.
Returns an `ODESolution`.
"""
function solve(
    solver::RK45Solver,
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

    # Dormand-Prince Butcher Coefficients
    c2 = 1.0 / 5.0
    c3 = 3.0 / 10.0
    c4 = 4.0 / 5.0
    c5 = 8.0 / 9.0
    c6 = 1.0
    c7 = 1.0

    a21 = 1.0 / 5.0
    a31 = 3.0 / 40.0;        a32 = 9.0 / 40.0
    a41 = 44.0 / 45.0;       a42 = -56.0 / 15.0;     a43 = 32.0 / 9.0
    a51 = 19372.0 / 6561.0;  a52 = -25360.0 / 2187.0; a53 = 64448.0 / 6561.0; a54 = -212.0 / 729.0
    a61 = 9017.0 / 3168.0;   a62 = -355.0 / 33.0;    a63 = 46732.0 / 5247.0;  a64 = 49.0 / 176.0;   a65 = -5103.0 / 18656.0
    a71 = 35.0 / 384.0;      a72 = 0.0;              a73 = 500.0 / 1113.0;    a74 = 125.0 / 192.0;  a75 = -2187.0 / 6784.0; a76 = 11.0 / 84.0

    # 5th-order weights
    b1 = 35.0 / 384.0
    b3 = 500.0 / 1113.0
    b4 = 125.0 / 192.0
    b5 = -2187.0 / 6784.0
    b6 = 11.0 / 84.0

    # Error weights (b_i - b̂_i)
    e1 = 71.0 / 57600.0
    e3 = -71.0 / 16695.0
    e4 = 71.0 / 1920.0
    e5 = -17253.0 / 339200.0
    e6 = 22.0 / 525.0
    e7 = -1.0 / 40.0

    k1 = f(u, p, t)
    n_fevals = 1
    n_steps = 0
    n_rejected = 0

    save_idx = 1
    has_save_times = save_times !== nothing && length(save_times) > 0
    if has_save_times && save_times[1] == t0
        save_idx = 2
    end

    while t < tf && n_steps < solver.max_steps
        if t + h > tf
            h = tf - t
        end

        # Stage 2
        u2 = u .+ h .* (a21 .* k1)
        k2 = f(u2, p, t + c2 * h)

        # Stage 3
        u3 = u .+ h .* (a31 .* k1 .+ a32 .* k2)
        k3 = f(u3, p, t + c3 * h)

        # Stage 4
        u4 = u .+ h .* (a41 .* k1 .+ a42 .* k2 .+ a43 .* k3)
        k4 = f(u4, p, t + c4 * h)

        # Stage 5
        u5 = u .+ h .* (a51 .* k1 .+ a52 .* k2 .+ a53 .* k3 .+ a54 .* k4)
        k5 = f(u5, p, t + c5 * h)

        # Stage 6
        u6 = u .+ h .* (a61 .* k1 .+ a62 .* k2 .+ a63 .* k3 .+ a64 .* k4 .+ a65 .* k5)
        k6 = f(u6, p, t + c6 * h)

        # 5th-order candidate solution
        u_next = u .+ h .* (b1 .* k1 .+ b3 .* k3 .+ b4 .* k4 .+ b5 .* k5 .+ b6 .* k6)

        # Stage 7 (FSAL evaluation)
        k7 = f(u_next, p, t + h)
        n_fevals += 6

        # Local error vector
        err_vec = h .* (e1 .* k1 .+ e3 .* k3 .+ e4 .* k4 .+ e5 .* k5 .+ e6 .* k6 .+ e7 .* k7)

        # Scaled error norm (RMS)
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
            k1_step = copy(k1)
            h_step = h
            t += h
            u = u_next
            k1 = k7 # FSAL reuse
            n_steps += 1

            if has_save_times
                while save_idx <= length(save_times) && save_times[save_idx] <= t + 1e-12
                    target_t = save_times[save_idx]
                    theta = clamp((target_t - t_old) / h_step, 0.0, 1.0)
                    # High-accuracy cubic Hermite polynomial interpolation (O(h^4))
                    h00 = 2.0 * theta^3 - 3.0 * theta^2 + 1.0
                    h10 = theta^3 - 2.0 * theta^2 + theta
                    h01 = -2.0 * theta^3 + 3.0 * theta^2
                    h11 = theta^3 - theta^2
                    u_interp = h00 .* u_old .+ h01 .* u .+ h_step .* (h10 .* k1_step .+ h11 .* k7)
                    push!(t_history, target_t)
                    push!(u_history, u_interp)
                    save_idx += 1
                end
            else
                push!(t_history, t)
                push!(u_history, copy(u))
            end

            # Step size adaptation for accepted step
            fac = solver.safety_factor * (1.0 / max(err_norm, 1e-5))^0.2
            fac = clamp(fac, 0.2, 5.0)
            h = clamp(h * fac, solver.min_step, solver.max_step)
        else
            # Step rejected
            n_rejected += 1
            fac = solver.safety_factor * (1.0 / err_norm)^0.25
            fac = clamp(fac, 0.1, 0.5)
            h = max(h * fac, solver.min_step)
        end
    end

    retcode = (t >= tf - 1e-10) ? :Success : :MaxStepsExceeded
    return ODESolution(t_history, u_history, retcode, n_steps, n_fevals, n_rejected)
end
