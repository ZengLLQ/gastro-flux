"""
    NCA: Non-Compartmental Pharmacokinetic Analysis Suite
"""

"""
    NCAMetrics

Standard non-compartmental pharmacokinetic parameters:
- `cmax`: Peak plasma concentration (mg/L)
- `tmax`: Time to peak concentration (h)
- `auc_last`: Area under curve to last observed point (mg*h/L)
- `auc_inf`: Area under curve extrapolated to infinity (mg*h/L)
- `aumc_inf`: Area under first moment curve extrapolated to infinity (mg*h²/L)
- `mrt`: Mean residence time (h) = AUMC / AUC
- `cl_f`: Apparent oral clearance (L/h) = Dose / AUC_inf
- `vz_f`: Apparent terminal volume of distribution (L) = (CL/F) / kel
- `t_half`: Terminal elimination half-life (h) = ln(2) / kel
"""
struct NCAMetrics
    cmax::Float64
    tmax::Float64
    auc_last::Float64
    auc_inf::Float64
    aumc_inf::Float64
    mrt::Float64
    cl_f::Float64
    vz_f::Float64
    t_half::Float64
end

"""
    compute_nca(t::Vector{Float64}, c::Vector{Float64}, dose_mg::Float64)

Performs robust trapezoidal non-compartmental analysis on a plasma concentration-time profile.
"""
function compute_nca(t::Vector{Float64}, c::Vector{Float64}, dose_mg::Float64)
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have equal length"
    @assert n >= 3 "At least 3 points required for NCA"

    cmax, idx = findmax(c)
    tmax = t[idx]

    auc_last = 0.0
    aumc_last = 0.0
    for i in 1:(n - 1)
        dt = t[i + 1] - t[i]
        auc_last += 0.5 * (c[i] + c[i + 1]) * dt
        aumc_last += 0.5 * (t[i] * c[i] + t[i + 1] * c[i + 1]) * dt
    end

    # Tail log-linear regression
    start_tail = max(1, floor(Int, 0.70 * n))
    tail_t = t[start_tail:end]
    tail_c = c[start_tail:end]

    kel = 0.1
    valid = tail_c .> 1e-6
    if count(valid) >= 3
        log_c = log.(tail_c[valid])
        sub_t = tail_t[valid]
        m_t = mean(sub_t)
        m_y = mean(log_c)
        denom = sum((sub_t .- m_t) .^ 2)
        if denom > 1e-12
            slope = sum((sub_t .- m_t) .* (log_c .- m_y)) / denom
            if slope < -1e-6
                kel = -slope
            end
        end
    end

    clast = c[end]
    tlast = t[end]

    auc_inf = auc_last + clast / kel
    aumc_inf = aumc_last + (clast * tlast / kel) + (clast / (kel^2))

    mrt = aumc_inf / max(auc_inf, 1e-12)
    cl_f = dose_mg / max(auc_inf, 1e-12)
    vz_f = cl_f / max(kel, 1e-12)
    t_half = log(2.0) / kel

    return NCAMetrics(cmax, tmax, auc_last, auc_inf, aumc_inf, mrt, cl_f, vz_f, t_half)
end
