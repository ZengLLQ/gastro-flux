using Test
using LinearAlgebra
using GastroFlux

@testset "High-Accuracy ODE Solvers" begin
    # 1. Exponential Decay: du/dt = -k*u => u(t) = u0 * exp(-k*t)
    k_decay = 2.5
    u0 = [10.0]
    tspan = (0.0, 3.0)
    f_decay(u, p, t) = -k_decay .* u

    solver_rk45 = RK45Solver(reltol=1e-8, abstol=1e-10)
    sol_rk45 = solve(solver_rk45, f_decay, u0, tspan)

    @test sol_rk45.retcode == :Success
    @test length(sol_rk45.t) > 5

    # Check accuracy at final time
    u_exact_final = u0[1] * exp(-k_decay * tspan[2])
    u_num_final = sol_rk45.u[end][1]
    rel_error = abs(u_num_final - u_exact_final) / u_exact_final
    @test rel_error < 1e-6

    # 2. Harmonic Oscillator: [x, v]' = [v, -omega^2 * x]
    omega = 3.0
    u0_ho = [1.0, 0.0] # x(0)=1, v(0)=0 => x(t) = cos(omega*t), v(t) = -omega*sin(omega*t)
    f_ho(u, p, t) = [u[2], -omega^2 * u[1]]

    sol_ho = solve(solver_rk45, f_ho, u0_ho, (0.0, 2.0 * π / omega))
    @test sol_ho.retcode == :Success
    # At t = 2*pi/omega, system returns to initial state
    @test isapprox(sol_ho.u[end][1], 1.0, atol=1e-4)
    @test isapprox(sol_ho.u[end][2], 0.0, atol=1e-4)

    # 3. Rosenbrock (ROS2) Stiff Solver: du/dt = -500.0 * (u - cos(t)) - sin(t)
    # Analytical solution: u(t) = cos(t) + (u0 - 1) * exp(-500*t)
    lambda_stiff = 500.0
    f_stiff(u, p, t) = [-lambda_stiff * (u[1] - cos(t)) - sin(t)]
    u0_stiff = [0.0] # u(0) = 0 => fast initial transient decaying to cos(t)

    solver_ros = RosenbrockSolver(reltol=1e-4, abstol=1e-6)
    sol_ros = solve(solver_ros, f_stiff, u0_stiff, (0.0, 2.0))

    @test sol_ros.retcode == :Success
    # At t = 2.0, exp(-500*2) is zero, so u(2) = cos(2.0)
    u_exact_stiff = cos(2.0)
    @test isapprox(sol_ros.u[end][1], u_exact_stiff, atol=2e-3)

    # 4. Dense interpolation / save_times matching
    target_times = [0.0, 0.5, 1.0, 1.5, 2.0]
    sol_dense = solve(solver_rk45, f_decay, u0, (0.0, 2.0); save_times=target_times)
    @test length(sol_dense.t) == length(target_times)
    @test all(isapprox.(sol_dense.t, target_times, atol=1e-10))
    for i in 1:length(target_times)
        u_exp = u0[1] * exp(-k_decay * target_times[i])
        @test isapprox(sol_dense.u[i][1], u_exp, rtol=1e-4)
    end
end
