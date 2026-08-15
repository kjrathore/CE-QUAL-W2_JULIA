using Test
using W2J
using Plots
using LinearAlgebra
using Random

"""
Literal transcription of `SUBROUTINE TRIDIAG` in `transport.f90:572-593`
(the real Fortran source, not re-derived from memory) -- kept unnormalized
until the final division exactly like the original, deliberately NOT
algebraically simplified to match `thomas_solve!`'s incrementally-normalized
form. Comparing against this, not just against `thomas_solve!` agreeing with
itself, is the actual validation: it proves the two different-looking
recurrences compute the same answer, which is the property that matters for
trusting `thomas_solve!` as a stand-in for `TRIDIAG`.
"""
function tridiag_fortran_literal(a, v, c, d)
    n = length(v)
    bta = similar(v, Float64)
    gma = similar(d, Float64)
    u = similar(d, Float64)
    bta[1] = v[1]
    gma[1] = d[1]
    for i in 2:n
        bta[i] = v[i] - a[i] / bta[i-1] * c[i-1]
        gma[i] = d[i] - a[i] / bta[i-1] * gma[i-1]
    end
    u[n] = gma[n] / bta[n]
    for i in (n-1):-1:1
        u[i] = (gma[i] - c[i] * u[i+1]) / bta[i]
    end
    return u
end

"Dense-matrix cross-check: builds the full tridiagonal matrix and solves with Julia's `\\`."
function dense_solve(a, b, c, d)
    n = length(b)
    A = zeros(Float64, n, n)
    for i in 1:n
        A[i, i] = b[i]
        i > 1 && (A[i, i-1] = a[i])
        i < n && (A[i, i+1] = c[i])
    end
    return A \ d
end

@testset "W2J" begin
    @testset "Hydrodynamics/Density: density" begin
        # Known physical reference points for fresh water (not just internal
        # self-consistency) -- standard textbook density-vs-temperature values.
        @test W2J.density(0.0, 0.0, 0.0, true, false, false) ≈ 999.8426 atol=1e-3
        @test W2J.density(4.0, 0.0, 0.0, true, false, false) ≈ 999.9750 atol=1e-3
        @test W2J.density(20.0, 0.0, 0.0, true, false, false) ≈ 998.2063 atol=1e-3
        @test W2J.density(25.0, 0.0, 0.0, true, false, false) ≈ 997.0480 atol=1e-3

        # Density must peak near 4 degC (fresh water's known density maximum) --
        # checks the polynomial shape, not just point values.
        rho_3 = W2J.density(3.0, 0.0, 0.0, true, false, false)
        rho_4 = W2J.density(4.0, 0.0, 0.0, true, false, false)
        rho_5 = W2J.density(5.0, 0.0, 0.0, true, false, false)
        @test rho_4 > rho_3 && rho_4 > rho_5

        # fresh_water/salt_water/susp_solids gate their respective terms --
        # turning them off must change the result relative to turning them on
        # (guards against the flags being silently ignored).
        @test W2J.density(10.0, 500.0, 0.0, true, false, false) != W2J.density(10.0, 500.0, 0.0, false, false, false)
        @test W2J.density(10.0, 500.0, 0.0, false, true, false) != W2J.density(10.0, 500.0, 0.0, false, false, false)
        @test W2J.density(10.0, 0.0, 100.0, false, false, true) != W2J.density(10.0, 0.0, 100.0, false, false, false)

        # T<=0 branch is a genuinely different polynomial (density.f90's IF/ELSE) --
        # exercise it explicitly, not just infer from the T>0 tests above.
        @test W2J.density(-1.0, 0.0, 0.0, true, false, false) ≈ 999.842594 atol=1e-6
    end

    @testset "Solvers/Tridiagonal: thomas_solve!" begin
        # Fixed 5x5 case (discrete 1D Laplacian-shaped system), cross-checked
        # against Julia's dense solver -- an independent, mature reference.
        a = [0.0, -1.0, -1.0, -1.0, -1.0]
        b = [2.0, 2.0, 2.0, 2.0, 2.0]
        c = [-1.0, -1.0, -1.0, -1.0, 0.0]
        d = [1.0, 0.0, 0.0, 0.0, 1.0]
        x = similar(d)
        W2J.thomas_solve!(a, b, c, d, x)
        @test x ≈ dense_solve(a, b, c, d)
        @test x ≈ [1.0, 1.0, 1.0, 1.0, 1.0]   # hand-verified: A*[1,1,1,1,1] == [1,0,0,0,1] == d

        # Randomized, diagonally-dominant systems (guarantees a well-conditioned,
        # numerically stable problem) of varying size, cross-checked against BOTH
        # the dense solver and a literal transcription of the real Fortran
        # TRIDIAG subroutine (transport.f90:572-593) -- see tridiag_fortran_literal
        # above. This is the actual "validate against TRIDIAG" the TODO asked for.
        rng = Random.Xoshiro(20260811)
        for trial in 1:200
            n = rand(rng, 2:50)
            b = 4.0 .+ rand(rng, n)                      # diagonal, dominant
            a = [0.0; -(1.0 .+ rand(rng, n - 1))]         # sub-diagonal, a[1] unused
            c = [-(1.0 .+ rand(rng, n - 1)); 0.0]         # super-diagonal, c[n] unused
            d = randn(rng, n)

            x = similar(d)
            W2J.thomas_solve!(a, b, c, d, x)

            @test x ≈ dense_solve(a, b, c, d) rtol=1e-9
            @test x ≈ tridiag_fortran_literal(a, b, c, d) rtol=1e-9
        end
    end

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
        @test geom.T2I == [6.7] && geom.ICEI == [0.0]
        @test geom.WTYPEC == ["FRESH"] && geom.GRIDC == ["RECT"]

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

    @testset "Core/InitGeometry: init_geometry!" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.InputReader.allocate_geometry!(g, geom)
        W2J.BathymetryReader.read_bathymetry!(geom, g, joinpath(detroit, "InputFiles", "bth1.csv"), 1; debug=false)

        W2J.init_geometry!(g, geom)

        # Bottom of the grid must equal ELBOT at every real segment (init-geom.F90:21).
        @test all(geom.EL[g.KMX, i] ≈ geom.ELBOT[1] for i in 2:30)

        # KB (bottom active layer) matches the documented range from CLAUDE.md
        # ("Active layer counts per branch: BR1 segments 2-11 (63-110 layers)").
        @test all(63 <= g.KB[i] <= 111 for i in 2:11)

        # KTI is identical for every segment -- true whenever ZERO_SLOPE(JW),
        # since every segment then shares the exact same EL(:,i) profile.
        @test g.ZERO_SLOPE == [true]
        @test length(unique(g.KTI[2:30])) == 1

        # No NaN/Inf and no negative volume anywhere in the real (non-padding) grid.
        @test all(isfinite, g.VOL[:, 2:30])
        @test all(>=(0.0), g.VOL[:, 2:30])
        @test all(isfinite, geom.DEPTHB[:, 2:30])

        # geom.H2 must carry the baseline H(:,JW) copy (input.F90:2265,2275) into
        # every layer, not just KT -- guards the class of bug where VOL/H2 look
        # "finite and non-negative" (0.0 satisfies both) while actually being
        # silently zero below the top active layer. Check a real wet segment
        # below its top layer has both nonzero H2 and nonzero VOL.
        wet_seg = g.CUS[1]
        below_kt = g.KTWB[1] + 1
        @test geom.H2[below_kt, wet_seg] > 0.0
        @test g.VOL[below_kt, wet_seg] > 0.0
        @test sum(g.VOL[:, wet_seg]) > g.VOL[g.KTWB[1], wet_seg]

        # Branch 1's upstream segments are shallower than the Jan-1 drawdown pool
        # (KB < KTI, i.e. dry) -- the model's own CUS logic must detect this and
        # mark the first properly-active segment accordingly (init-geom.F90:474-484).
        @test g.CUS[1] > g.US[1]
        @test all(g.KB[i] < g.KTI[i] for i in g.US[1]:(g.CUS[1]-1))
        @test all(g.KB[i] >= g.KTI[i] for i in g.CUS[1]:g.DS[1])

        # JBUH/JBDH (branch-to-branch linkage, computed via Core/Grid.jl's
        # BranchNetwork inside compute_bottom_layers!) match Detroit's known
        # topology: branches 2-4 connect downstream into branch 1.
        @test geom.JBDH[2] == 1 && geom.JBDH[3] == 1 && geom.JBDH[4] == 1

        # TRAPEZOIDAL(JW) computed unconditionally by init_geometry! (only
        # needs GRIDC, no Tier-1 dependency) -- Detroit is GRIDC="RECT".
        @test geom.TRAPEZOIDAL == [false]

        # FRESH_WATER/SALT_WATER NOT computed automatically by init_geometry!
        # (needs constituents/tds_active, not read by InputReader.jl yet) --
        # must stay empty until explicitly requested.
        @test isempty(geom.FRESH_WATER) && isempty(geom.SALT_WATER)

        # compute_water_type_flags! respects constituents/tds_active as real
        # gates, not decoration -- turning either off must suppress
        # FRESH_WATER even though WTYPEC=="FRESH" (init.F90:172-173).
        W2J.compute_water_type_flags!(g, geom; constituents=true, tds_active=true)
        @test geom.FRESH_WATER == [true] && geom.SALT_WATER == [false]
        W2J.compute_water_type_flags!(g, geom; constituents=false, tds_active=true)
        @test geom.FRESH_WATER == [false]
        W2J.compute_water_type_flags!(g, geom; constituents=true, tds_active=false)
        @test geom.FRESH_WATER == [false]

        # compute_areas_volumes! must refuse (not silently compute wrong RECT-
        # formula areas) for a waterbody flagged TRAPEZOIDAL, since that branch
        # isn't ported -- flip Detroit's own (real) geom.TRAPEZOIDAL to prove
        # the guard actually fires, not just that it exists unreachable.
        geom.TRAPEZOIDAL = [true]
        @test_throws ErrorException W2J.compute_areas_volumes!(g, geom)
    end

    @testset "Hydrodynamics/FreeSurface: hydrodynamic_step! (zero-flow sanity check)" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.InputReader.allocate_geometry!(g, geom)
        W2J.BathymetryReader.read_bathymetry!(geom, g, joinpath(detroit, "InputFiles", "bth1.csv"), 1; debug=false)
        g, geom, net = W2J.init_geometry!(g, geom)
        W2J.compute_dlxrho!(g, geom)
        W2J.allocate_hydro_state!(g)

        elws_before = copy(geom.ELWS)
        dlt = tc.DLTMAX[1]

        for step in 1:5
            W2J.hydrodynamic_step!(g, geom, net, dlt)
        end

        # Zero boundary flow, uniform density (single T2I), zero slope (SINAC=0
        # for Detroit) -- every forcing term in the reduced-physics solve is
        # genuinely zero, so water surface elevation must stay at its initial
        # value (not just "finite") and velocity must stay at 0. A bug in branch
        # ordering, the BHRHO/A/V/C tridiagonal coefficients, or DH_INTERNAL
        # boundary coupling would show up as drift or NaN here, not as a
        # trivially-passing check.
        @test all(isfinite, geom.ELWS[2:30])
        @test all(isfinite, g.U[:, 2:30])
        @test maximum(abs, geom.ELWS[2:30] .- elws_before[2:30]) < 1e-9
        @test maximum(abs, g.U[:, 2:30]) < 1e-9
    end

    @testset "Hydrodynamics/FreeSurface: long-run stability (2000 steps)" begin
        # The 5-step check above proves the tridiagonal assembly/branch
        # sequencing is right at start-up; it does NOT prove the reduced-
        # physics solve stays stable over a real simulation length. Detroit's
        # own control file spans a full year (TMSTRT=1, TMEND=365) at
        # DLTMAX=1200s -- ~26,208 steps. 2000 steps here (not the full
        # 26,208) keeps the test suite fast while still covering two orders
        # of magnitude more steps than the smoke-test-sized check above; a
        # slow floating-point drift that only shows up after hundreds of
        # steps (not 5) would be caught here.
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.InputReader.allocate_geometry!(g, geom)
        W2J.BathymetryReader.read_bathymetry!(geom, g, joinpath(detroit, "InputFiles", "bth1.csv"), 1; debug=false)
        g, geom, net = W2J.init_geometry!(g, geom)
        W2J.compute_dlxrho!(g, geom)
        W2J.allocate_hydro_state!(g)

        elws_before = copy(geom.ELWS)
        dlt = tc.DLTMAX[1]
        for _ in 1:2000
            W2J.hydrodynamic_step!(g, geom, net, dlt)
        end

        @test all(isfinite, geom.ELWS[2:30])
        @test all(isfinite, g.U[:, 2:30])
        @test maximum(abs, geom.ELWS[2:30] .- elws_before[2:30]) < 1e-8
        @test maximum(abs, g.U[:, 2:30]) < 1e-8
    end

    @testset "Core/Parallel: parallel_foreach" begin
        # threshold=0 forces the Threads.@threads branch; threshold=typemax
        # forces the serial branch -- both must visit every index exactly
        # once and agree with each other, regardless of how many threads
        # the test process actually has (correctness of the dispatch logic
        # itself, not a timing test).
        for range in (1:1, 1:10, 1:200)
            hit_threaded = zeros(Int, length(range))
            W2J.parallel_foreach(range; threshold=0) do i
                hit_threaded[i - first(range) + 1] += 1
            end
            @test all(==(1), hit_threaded)

            hit_serial = zeros(Int, length(range))
            W2J.parallel_foreach(range; threshold=typemax(Int)) do i
                hit_serial[i - first(range) + 1] += 1
            end
            @test hit_serial == hit_threaded
        end
    end

    @testset "Simulation: run_zero_flow_sanity_check! end-to-end TSR output" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        outdir = joinpath(detroit, "Outputs")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.InputReader.allocate_geometry!(g, geom)
        W2J.BathymetryReader.read_bathymetry!(geom, g, joinpath(detroit, "InputFiles", "bth1.csv"), 1; debug=false)
        g, geom, net = W2J.init_geometry!(g, geom)
        W2J.compute_dlxrho!(g, geom)
        W2J.allocate_hydro_state!(g)

        # One representative wet segment per branch (branch 1's first properly-
        # active segment CUS[1], plus the downstream-most segment of the other
        # three branches) -- not a hardcoded Detroit assumption baked into
        # Simulation.jl itself, just this test's choice of what to sample.
        output_segments = [g.CUS[1], g.DS[2], g.DS[3], g.DS[4]]

        g, geom, jday_final = W2J.run_zero_flow_sanity_check!(g, geom, net, tc;
            nsteps=5, output_dir=outdir, output_segments=output_segments, base_name="tsr_zeroflow")

        @test jday_final ≈ tc.TMSTRT + 5 * tc.DLTMAX[1] / 86400.0

        for seg in output_segments
            path = joinpath(outdir, "tsr_zeroflow_seg$(seg).csv")
            @test isfile(path)
            lines = readlines(path)
            @test lines[1] == "JDAY,DLT(s),ELWS(m),U(ms-1)"
            @test length(lines) == 1 + 6   # header + initial condition + 5 steps

            # Zero-flow sanity check: ELWS must not drift, U must stay at 0.
            elws_vals = [parse(Float64, split(l, ',')[3]) for l in lines[2:end]]
            u_vals = [parse(Float64, split(l, ',')[4]) for l in lines[2:end]]
            @test maximum(abs, elws_vals .- elws_vals[1]) < 1e-9
            @test all(iszero, u_vals)
        end
    end

    @testset "Core/Grid: build_branch_network" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)

        net = W2J.build_branch_network(g)

        # Detroit topology (CLAUDE.md validated values): US=[2,14,22,28],
        # DS=[11,19,25,30], UHS=[0,0,0,0] (all external), DHS=[0,9,10,11]
        # (branches 2-4 connect into branch 1 at those segments).
        @test net.segment_range == [2:11, 14:19, 22:25, 28:30]
        @test net.waterbody_of_branch == [1, 1, 1, 1]
        @test net.upstream_branch == [0, 0, 0, 0]           # UHS all 0 -- no internal upstream connections
        @test net.downstream_branch == [0, 1, 1, 1]          # branch 1 -> external (dam); 2,3,4 -> branch 1
        @test net.downstream_waterbody == [0, 1, 1, 1]
        @test net.segment_to_branch[9] == 1 && net.segment_to_branch[10] == 1 && net.segment_to_branch[11] == 1
        @test net.segment_to_branch[15] == 2   # a segment inside branch 2's own range
    end

    @testset "Core/Grid: branch_processing_order" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.compute_boundary_flags!(g, geom)
        net = W2J.build_branch_network(g)

        # Detroit: branches 2-4 each depend on branch 1 (DH_INTERNAL). Detroit's
        # own branch numbering already happens to be topological, so this alone
        # doesn't prove generality -- see the reversed-topology case below.
        @test W2J.branch_processing_order(g, net, 1) == [1, 2, 3, 4]

        # Synthetic REVERSED topology: branch 4 is the terminal (downstream-
        # most, externally-draining) branch; branches 1-3 all connect DOWN
        # into it. Numeric order 1,2,3,4 would be WRONG here -- a correct
        # general (not hardcoded-to-ascending-numbering) algorithm must still
        # place 4 before 1, 2, and 3.
        g2 = W2J.W2Core.W2Global()
        g2.NBR = 4; g2.NWB = 1; g2.IMX = 42
        g2.US = [2, 12, 22, 32]; g2.DS = [10, 20, 30, 40]
        g2.BS = [1]; g2.BE = [4]
        g2.UHS = [0, 0, 0, 0]      # all external upstream
        g2.DHS = [32, 32, 32, 0]  # branches 1-3 -> branch 4 (segment 32); branch 4 -> external
        g2.JBDN = [4]
        geom2 = W2J.W2Core.W2Geometry()
        geom2.SLOPE = zeros(4)   # compute_boundary_flags! needs SLOPE for ZERO_SLOPE, unrelated to this test
        W2J.compute_boundary_flags!(g2, geom2)
        net2 = W2J.build_branch_network(g2)

        @test net2.downstream_branch == [4, 4, 4, 0]
        order2 = W2J.branch_processing_order(g2, net2, 1)
        @test findfirst(==(4), order2) < findfirst(==(1), order2)
        @test findfirst(==(4), order2) < findfirst(==(2), order2)
        @test findfirst(==(4), order2) < findfirst(==(3), order2)
    end

    @testset "Core/Grid: branch_processing_tiers (parallel-processing groundwork)" begin
        detroit = joinpath(@__DIR__, "..", "..", "DetroitReservoir")
        g, geom, tc = W2J.InputReader.read_control_file(joinpath(detroit, "w2_con.csv"); debug=false)
        W2J.compute_boundary_flags!(g, geom)
        net = W2J.build_branch_network(g)

        # Branch 1 has no internal dependency -> tier 1 alone. Branches 2-4 each
        # depend only on branch 1 (already resolved) -> tier 2 together, the
        # concrete case this function exists for: three branches genuinely safe
        # to process under Threads.@threads concurrently.
        tiers = W2J.branch_processing_tiers(g, net, 1)
        @test tiers == [[1], [2, 3, 4]]

        # Flattening tiers must reproduce branch_processing_order exactly.
        @test reduce(vcat, tiers) == W2J.branch_processing_order(g, net, 1)

        # Synthetic reversed topology (same construction as the order test
        # above): branch 4 is the true downstream/terminal branch; 1-3 all
        # depend only on it -> tier 1 is [4] alone, tier 2 is [1,2,3] together.
        g2 = W2J.W2Core.W2Global()
        g2.NBR = 4; g2.NWB = 1; g2.IMX = 42
        g2.US = [2, 12, 22, 32]; g2.DS = [10, 20, 30, 40]
        g2.BS = [1]; g2.BE = [4]
        g2.UHS = [0, 0, 0, 0]
        g2.DHS = [32, 32, 32, 0]
        g2.JBDN = [4]
        geom2 = W2J.W2Core.W2Geometry()
        geom2.SLOPE = zeros(4)
        W2J.compute_boundary_flags!(g2, geom2)
        net2 = W2J.build_branch_network(g2)

        tiers2 = W2J.branch_processing_tiers(g2, net2, 1)
        @test tiers2 == [[4], [1, 2, 3]]
    end
end
