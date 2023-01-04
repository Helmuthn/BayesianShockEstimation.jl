using BayesianShockEstimation

@testset "riemann_solver_burgers" begin
    @testset "zeros" begin
        ul = 0
        ur = 0
        flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
        @test flux == 0
    end

    @testset "ul > 0 and ur > 0" begin
        ul = 1
        ur = 2
        flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
        truth = BayesianShockEstimation.flux_burgers(ul)
        @test flux == truth
    end

    @testset "ul < 0 and ur < 0" begin
        ul = -1
        ur = -2
        flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
        truth = BayesianShockEstimation.flux_burgers(ur)
        @test flux == truth
    end

    @testset "ul < 0 and ur > 0" begin
        ul = -1
        ur = 2
        flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
        truth = 0
        @test flux == truth
    end

    @testset "ul > 0 and ur < 0" begin
        @testset "Zero Speed Shock" begin
            # Not well defined, either side is the same
            ul = 1
            ur = -1
            flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
            truth = BayesianShockEstimation.flux_burgers(ul)
            @test flux == truth
        end

        @testset "Left Moving Shock" begin
            ul = 1
            ur = -2
            flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
            truth = BayesianShockEstimation.flux_burgers(ur)
            @test flux == truth
        end

        @testset "Right Moving Shock" begin
            ul = 2
            ur = -1
            flux = BayesianShockEstimation.riemann_solver_burgers(ul, ur)    
            truth = BayesianShockEstimation.flux_burgers(ul)
            @test flux == truth
        end
    end
end
