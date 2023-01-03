"""
Sample from RTS smoother backwards sampler for random walk
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
Helper Function computing variance for RTS sampler

Returns forward estimates with RTS variances
"""
function RTS_Smooth(observations, σ_w, σ_v)
    estimates, variances = KalmanFilterWalk(observations, σ_w, σ_v)
    RTSvariances = BackwardsVarianceWalk(variances, σ_w)
    return estimates, variances, RTSvariances
end
export RTS_Smooth

"""
Helper function computing Kalman filter for the given model
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
Helper Function computing variance from RTS smoother
"""
function BackwardsVarianceWalk(variances, σ_w)
    RTSvariances = zeros(length(variances))
    for i in length(variances)-1:-1:1
        weight = σ_w^2/(variances[i] + σ_w^2)
        RTSvariances[i] = weight * variances[i] + (1-weight)^2 * RTSvariances[i+1]
    end
    return RTSvariances
end

