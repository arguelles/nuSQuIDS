using Test
using nuSQuIDS

@testset "nuSQuIDS.jl" begin
    @testset "Constants" begin
        units = Const()
        @test GeV(units) > 0
        @test km(units) > 0
        @test GeV(units) == 1e9 * eV(units)
        @test km(units) == 1000 * meter(units)
    end

    @testset "Vacuum Oscillations" begin
        units = Const()

        # Create single energy nuSQUIDS
        nus = nuSQUIDS(3, neutrino)

        # Set vacuum body and track
        Set_Body(nus, Vacuum())
        Set_Track(nus, VacuumTrack(1000.0 * km(units)))

        # Set energy
        Set_E(nus, 1.0 * GeV(units))

        # Set mixing parameters (standard values)
        Set_MixingAngle(nus, 0, 1, 0.5836)  # theta_12
        Set_MixingAngle(nus, 0, 2, 0.1496)  # theta_13
        Set_MixingAngle(nus, 1, 2, 0.8587)  # theta_23
        Set_SquareMassDifference(nus, 1, 7.5e-5 * eV(units)^2)  # dm21^2
        Set_SquareMassDifference(nus, 2, 2.5e-3 * eV(units)^2)  # dm31^2

        # Set initial state: pure nu_mu
        Set_initial_state(nus, [0.0, 1.0, 0.0], flavor)

        # Evolve
        EvolveState(nus)

        # Check probabilities sum to 1
        p_e = EvalFlavor(nus, 0)
        p_mu = EvalFlavor(nus, 1)
        p_tau = EvalFlavor(nus, 2)

        @test p_e >= 0 && p_e <= 1
        @test p_mu >= 0 && p_mu <= 1
        @test p_tau >= 0 && p_tau <= 1
        @test isapprox(p_e + p_mu + p_tau, 1.0, atol=1e-10)
    end

    @testset "Energy Arrays" begin
        E_lin = linspace(1e9, 1e12, 10)
        @test length(E_lin) == 10
        @test E_lin[1] ≈ 1e9
        @test E_lin[end] ≈ 1e12

        E_log = logspace(1e9, 1e12, 10)
        @test length(E_log) == 10
        @test E_log[1] ≈ 1e9
        @test E_log[end] ≈ 1e12
    end

    @testset "SU_vector" begin
        # Test basic SU_vector operations
        v = SU_vector(3)
        @test Dim(v) == 3
        @test Size(v) == 9  # 3^2 components
    end

    @testset "Bodies" begin
        # Test body creation
        vac = Vacuum()
        @test density(vac, VacuumTrack(1.0)) == 0.0
        @test ye(vac, VacuumTrack(1.0)) == 0.0

        # Constant density
        cd = ConstantDensity(3.0, 0.5)
        @test density(cd, ConstantDensityTrack(1.0)) ≈ 3.0
        @test ye(cd, ConstantDensityTrack(1.0)) ≈ 0.5
    end

    @testset "Resource Path" begin
        path = getResourcePath()
        @test isdir(path) || !isempty(path)
    end
end

println("All tests passed!")
