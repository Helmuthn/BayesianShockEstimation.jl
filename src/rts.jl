"""
    RTS_sample(estimates, variances, RTSvariances, σ_w)

Sample from RTS smoother backwards sampler for random walk model.

Uses the Kalman filter estimates and variances, as well as the RTS variances
generated by `RTS_Smooth`.

# Arguments
 - `estimates`    : Kalman filter estimates
 - `variances`    : Kalman filter variances
 - `RTSvariances` : RTS smoother variances
 - `σ_w`          : System noise standard deviation


"""
function RTS_sample(estimates, variances, RTSvariances, σ_w)
    # Start by generating the noise
    out = sqrt.(RTSvariances) .* randn(length(estimates))
    
    # Add in the means
    for i in length(out)-1:-1:1
        weight = σ_w^2/(σ_w^2 + variances[i])
        out[i] += weight * estimates[i] + (1 - weight) * out[i+1]
    end

    return out
end

"""
    RTS_Smooth(observations, σ_w, σ_v)

Helper Function computing variance for RTS sampler by computing
the Kalman filter estimates as well as the RTS smoother variances.

Applies to the random walk model

```math
\\begin{aligned}
x_{t+1} &= x_t + w_t\\\\
y_t     &= x_t + v_t\\\\
\\end{aligned}
```

where ``w_t`` and ``v_t`` are centered Gaussian random variables
with variances ``σ_w^2`` and ``σ_v^2``, respectively.

# Arguments
 - `observations` : Timeseries observation
 - `σ_w`          : System noise standard deviation
 - `σ_v`          : Observaton noise standard deviation

# Returns
    `estimates, variances, RTSvariances`

 - `estimates`    : Kalman filter estimates
 - `variances`    : Kalman filter variances
 - `RTSvariances` : RTS smoother variances
"""
function RTS_Smooth(observations, σ_w, σ_v)
    estimates, variances = KalmanFilterWalk(observations, σ_w, σ_v)
    RTSvariances = BackwardsVarianceWalk(variances, σ_w)
    return estimates, variances, RTSvariances
end
export RTS_Smooth

"""
    KalmanFilterWalk(observations, σ_w, σ_v)

Computes the Kalman filter for a random walk model.

```math
\\begin{aligned}
x_{t+1} &= x_t + w_t\\\\
y_t     &= x_t + v_t\\\\
\\end{aligned}
```

where ``w_t`` and ``v_t`` are centered Gaussian random variables
with variances ``σ_w^2`` and ``σ_v^2``, respectively.

# Arguments
 - `observations` : Timeseries observation
 - `σ_w`          : System noise standard deviation
 - `σ_v`          : Observaton noise standard deviation

# Returns
    `estimates, variances`

 - `estimates` : Kalman filter estimates
 - `variances` : Kalman filter variances
"""
function KalmanFilterWalk(observations, σ_w, σ_v)
    estimates = zeros(length(observations))
    variances = zeros(length(observations))

    # Initialize the Kalman filter
    estimates[1] = observations[0]
    variances[1] = σ_v^2

    for i in 2:length(observations)
        weight = σ_v^2/(variances[i-1] + σ_v^2 + σ_w^2)
        estimates[i] = weight * estimates[i-1] + (1 - weight) * observations[i]
        variances[i] = weight * (variances[i-1] + σ_w^2)
    end

    return estimates, variances
end

"""
    BackwardsVarianceWalk(variances, σ_w)

Given Kalman filter variances, computes the RTS smoother variances.

# Arguments
 - `variances` : Kalman filter variances
 - `σ_w`       : System noise power

# Returns

The variances for the RTS smoother
"""
function BackwardsVarianceWalk(variances, σ_w)
    RTSvariances = zeros(length(variances))
    for i in length(variances)-1:-1:1
        weight = σ_w^2/(variances[i] + σ_w^2)
        RTSvariances[i] = weight * variances[i] + (1-weight)^2 * RTSvariances[i+1]
    end
    return RTSvariances
end

