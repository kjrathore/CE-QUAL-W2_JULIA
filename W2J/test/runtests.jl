using Test
using W2J
using Plots

@testset "W2J" begin
    # TODO: Solvers/Tridiagonal.jl's thomas_solve! is the first thing worth
    # unit-testing in isolation — it's pure math with no W2-specific
    # dependencies. A simple known-solution tridiagonal system (e.g. from a
    # textbook example) would validate it before any Fortran comparison.
    #
    # The REAL validation target is end-to-end: a single unbranched
    # reservoir test case compared against existing Fortran reference
    # output. See test/reference_cases/ (currently empty) and README
    # "Base Module / Walking-Skeleton Plan".

    @testset "IO end-to-end: read_control_file -> allocate -> read_bathymetry!" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        con_file = joinpath(detroit, "w2_con.csv")
        bth_file = joinpath(detroit, "InputFiles", "bth1.csv")
        outdir = joinpath(detroit, "Outputs")
        mkpath(outdir)

        g, geom, tc = W2J.InputReader.read_control_file(con_file; debug=false)

        # Validated Detroit values -- see CLAUDE.md "Validated values for Detroit".
        @test g.NWB == 1 && g.NBR == 4 && g.IMX == 31 && g.KMX == 117
        @test g.US == [2, 14, 22, 28]
        @test g.DS == [11, 19, 25, 30]
        @test g.BS == [1] && g.BE == [4]
        @test geom.ELBOT == [364.0]

        W2J.InputReader.allocate_geometry!(g, geom)
        @test size(geom.B) == (g.KMX, g.IMX)
        @test length(geom.DLX) == g.IMX

        bathy = W2J.BathymetryReader.read_bathymetry!(geom, g, bth_file, 1; debug=false)
        @test geom.DLX[2] ≈ 844.9
        @test geom.ELWS[2] ≈ 441.09
        @test bathy.H[1] ≈ 1.95

        p = W2J.LongitudinalProfile.plot_longitudinal_profile(g, geom, Dict(1 => bathy))
        @test p !== nothing
        savefig(p, joinpath(outdir, "detroit_longitudinal_profile.png"))
    end
end
