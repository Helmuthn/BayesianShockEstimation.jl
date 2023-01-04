using BayesianShockEstimation

@testset "KalmanFilterWalk" begin
    observations = zeros(2)
    σ_v = 1
    σ_w = 1
    estimates, variances = KalmanFilterWalk(observations, σ_w, σ_v)

    @test estimates[1] == observations[1]
    @test variances[1] == σ_v^2
    @test estimates[2] == observations[2]
    @test variances[2] ≈ σ_v^2 * (σ_v^2 + σ_w^2)/(2 * σ_v^2 + σ_w^2)
end

@testset "BackwardsVarianceWalk" begin
    variances = [1.0, 1.0]
    σ_w = 1.0
    RTSVariances = BackwardsVarianceWalk(variances, σ_w)
    @test RTSVariances[end] == variances[end]
    @test RTSVariances[1] == 0.75
end

@testset "RTS_Smooth" begin
    observations = [1,1]
    σ_w = 1
    σ_v = 1
    estimates, variances, RTSvariances = RTS_Smooth(observations, σ_w, σ_v)
    @test estimates == [1,1]
end
