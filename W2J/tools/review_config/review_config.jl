# ==============================================================================
# tools/review_config/review_config.jl
#
# Ported CE-QUAL-W2 preprocessor validation checks (from the reference Fortran
# preprocessor `v455/preprocessor_422/preprocessor_4.22/pre_ivf_422.f90`,
# ~13,725 lines) -- runs a SUBSET of those checks against `w2_config.yaml`
# (produced by tools/xlsx_to_yaml/convert.jl).
#
# SCOPE (deliberately not exhaustive): that Fortran file has roughly 1,174
# distinct validation sites (`CALL ERRORS`/`CALL WARNINGS`, each paired with a
# `WRITE(ERR,...)`/`WRITE(WRN,...)` message), spread across nearly every
# input section -- structures, withdrawals, gates, pipes, meteorology,
# kinetics rate coefficients, output control, and more. Porting all of them
# is a multi-session effort. This first pass covers only the sections already
# reflected in InputReader.jl / xlsx_to_yaml's output: dimensions/NPROC/
# CLOSEC, branch & grid geometry (including JBDN), and CST constituent
# activation. A clean run of this script means "the checks below passed," NOT
# "the control file is fully valid" -- see W2J_README.md "Preprocessor checks
# ported" for the running list of what is and isn't covered yet.
#
# Every check cites its Fortran source line range (as of pre_ivf_422.f90 in
# this repo) so it can be verified or extended directly against the original,
# rather than re-derived from memory.
# ==============================================================================

using YAML

struct Finding
    severity::Symbol   # :error or :warning
    section::String
    message::String
end

"Normalizes a config value that may be a bare scalar (NWB=1 case) or a YAML list (NWB>1)."
as_vector(x) = x isa AbstractVector ? x : Any[x]

"""
Finds a top-level config section by key prefix. `exclude` filters out keys
containing any of the given substrings -- e.g. "CST -" alone would also match
"CST DERI - ..." and "CST FLUX - ...", so those are excluded when looking for
plain CST.
"""
function find_section(cfg, prefix::String; exclude::Vector{String}=String[])
    for (k, v) in cfg
        k isa String || continue
        startswith(k, prefix) || continue
        any(e -> occursin(e, k), exclude) && continue
        return v
    end
    return nothing
end

# ------------------------------------------------------------------------------
# Dimensions / NPROC / CLOSEC -- pre_ivf_422.f90 lines 3499-3512
# ------------------------------------------------------------------------------
function check_dimensions(cfg)
    findings = Finding[]
    dims = find_section(cfg, "GRID/NPROC")
    dims === nothing && return findings

    nproc = dims["NPROC"]
    if nproc < 0
        push!(findings, Finding(:error, "GRID/NPROC", "NPROC=$nproc < 0, must be > 0"))
    elseif nproc == 0
        push!(findings, Finding(:error, "GRID/NPROC", "NPROC=0, must be > 0 (W2 sets NPROC=1 when NPROC=0)"))
    elseif nproc > 4
        push!(findings, Finding(:warning, "GRID/NPROC", "NPROC=$nproc > 4 -- make sure this does not lead to long run times"))
    end

    closec = strip(string(dims["CLOSEC"]))
    if !(uppercase(closec) in ("OFF", "ON"))
        push!(findings, Finding(:error, "GRID/NPROC", "CLOSEC must be \"OFF\" or \"ON\", got \"$closec\""))
    end

    return findings
end

# ------------------------------------------------------------------------------
# Branch & grid geometry -- pre_ivf_422.f90 lines 2794-2820 (bathymetry/branch
# geometry), 3306-3335 (UHS/DHS validity), 3489-3494 (slope), 8417-8427 (JBDN)
# ------------------------------------------------------------------------------
function check_branch_geometry(cfg)
    findings = Finding[]
    grid = find_section(cfg, "BRANCH GRID")
    loc = find_section(cfg, "LOCATION")
    dims = find_section(cfg, "GRID/NPROC")
    (grid === nothing || loc === nothing || dims === nothing) && return findings

    US = as_vector(grid["US"]); DS = as_vector(grid["DS"])
    UHS = as_vector(grid["UHS"]); DHS = as_vector(grid["DHS"])
    SLOPE = as_vector(grid["SLOPE"]); SLOPEC = as_vector(grid["SLOPEC"])
    NBR = dims["NBR"]; IMX = dims["IMX"]
    BS = as_vector(loc["BS"]); BE = as_vector(loc["BE"]); JBDN = as_vector(loc["JBDN"])

    # pre_ivf_422.f90:2794-2802 -- sum(BE-BS+1) over waterbodies must equal NBR
    isum = sum(BE[jw] - BS[jw] + 1 for jw in eachindex(BS))
    if isum != NBR
        push!(findings, Finding(:error, "BRANCH GRID",
            "sum(BE-BS+1) over waterbodies = $isum, does not equal NBR=$NBR (BS/BE error)"))
    end

    # pre_ivf_422.f90:2806-2809 -- US(1) must equal 2 (first segment is inactive)
    if US[1] != 2
        push!(findings, Finding(:error, "BRANCH GRID", "US(1)=$(US[1]) must equal 2 -- first segment must be inactive"))
    end

    # pre_ivf_422.f90:2810-2816 -- exactly 2 inactive boundary segments between branches
    for jb in 1:(NBR-1)
        diff = US[jb+1] - DS[jb]
        if diff != 3
            push!(findings, Finding(:error, "BRANCH GRID",
                "US(branch $(jb+1))-DS(branch $jb) = $diff, expected 3 -- missing or extra inactive boundary segments"))
        end
    end

    # pre_ivf_422.f90:2817-2820 -- DS(last branch)+1 must equal IMX
    if DS[BE[end]] + 1 != IMX
        push!(findings, Finding(:error, "BRANCH GRID",
            "DS(branch $(BE[end]))+1 = $(DS[BE[end]]+1), does not equal IMX=$IMX"))
    end

    # pre_ivf_422.f90:3489-3494 -- corrected slope cannot exceed actual slope
    for jb in eachindex(SLOPE)
        if SLOPE[jb] < SLOPEC[jb]
            push!(findings, Finding(:error, "BRANCH GRID",
                "SLOPE(branch $jb)=$(SLOPE[jb]) < SLOPEC(branch $jb)=$(SLOPEC[jb]) -- corrected slope cannot exceed actual slope"))
        end
    end

    # pre_ivf_422.f90:3306-3335 -- DHS/UHS, if positive, must fall within another branch's US..DS range
    for jb in eachindex(DHS)
        dhs = DHS[jb]
        if dhs > 0 && !any(jjb -> jjb != jb && dhs >= US[jjb] && dhs <= DS[jjb], eachindex(US))
            push!(findings, Finding(:error, "BRANCH GRID", "DHS(branch $jb)=$dhs does not fall within any other branch's US..DS range"))
        end
    end
    for jb in eachindex(UHS)
        uhs = UHS[jb]
        if uhs > 0 && !any(jjb -> jjb != jb && uhs >= US[jjb] && uhs <= DS[jjb], eachindex(US))
            push!(findings, Finding(:error, "BRANCH GRID", "UHS(branch $jb)=$uhs does not fall within any other branch's US..DS range"))
        end
    end

    # pre_ivf_422.f90:8417-8427 -- JBDN(JW) must be one of waterbody JW's own branches
    for jw in eachindex(JBDN)
        if !(BS[jw] <= JBDN[jw] <= BE[jw])
            push!(findings, Finding(:error, "LOCATION",
                "JBDN(waterbody $jw)=$(JBDN[jw]) is not within its own branch range [$(BS[jw]),$(BE[jw])]"))
        end
    end

    return findings
end

# ------------------------------------------------------------------------------
# CST constituent activation -- pre_ivf_422.f90 lines 10383-10443
# ------------------------------------------------------------------------------
function check_cst(cfg)
    findings = Finding[]
    cst = find_section(cfg, "CST -"; exclude=["DERI", "FLUX"])
    cst === nothing && return findings

    for (name, c) in cst
        active = strip(uppercase(string(get(c, "active", ""))))
        if !(active in ("ON", "OFF"))
            push!(findings, Finding(:error, "CST", "constituent $name: active=\"$active\" must be \"ON\" or \"OFF\""))
        end

        # pre_ivf_422.f90:10409-10421 -- initial conc < -2 is an error; ==0 while
        # active is a warning, except for age/residence-time constituents. The
        # Fortran source checks its own CNAME1 array (a fixed display name, not
        # user-editable) for "Residence time" / "AGE" -- we don't carry that
        # third name variant, so check both the short name (CNAME2, e.g. "Gen1")
        # AND the long name (CNAME, e.g. "Water age, days") for "AGE"/"RESIDENCE".
        initial_conc = as_vector(get(c, "initial_conc", []))
        long_name_upper = uppercase(string(get(c, "long_name", "")))
        exempt = occursin("AGE", uppercase(string(name))) || occursin("AGE", long_name_upper) ||
                 occursin("RESIDENCE", long_name_upper)
        for (jw, v) in enumerate(initial_conc)
            v isa Number || continue
            if v < -2.0
                push!(findings, Finding(:error, "CST", "constituent $name: initial_conc(waterbody $jw)=$v < -2.0"))
            elseif active == "ON" && v == 0.0 && !exempt
                push!(findings, Finding(:warning, "CST", "constituent $name: initial_conc(waterbody $jw)=0 while active=ON"))
            end
        end

        # pre_ivf_422.f90:10428-10443 -- print control must be ON/OFF; printing
        # an inactive constituent is a warning, not an error.
        print_snp = as_vector(get(c, "print_snp", []))
        for (jw, v) in enumerate(print_snp)
            pv = strip(uppercase(string(v)))
            if !(pv in ("ON", "OFF"))
                push!(findings, Finding(:error, "CST", "constituent $name: print_snp(waterbody $jw)=\"$pv\" must be \"ON\" or \"OFF\""))
            elseif pv == "ON" && active != "ON"
                push!(findings, Finding(:warning, "CST", "constituent $name: print_snp(waterbody $jw)=ON but active=$active"))
            end
        end
    end

    return findings
end

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
function run_checks(cfg)
    findings = Finding[]
    append!(findings, check_dimensions(cfg))
    append!(findings, check_branch_geometry(cfg))
    append!(findings, check_cst(cfg))
    return findings
end

function main(yaml_path::AbstractString)
    cfg = YAML.load_file(yaml_path)
    findings = run_checks(cfg)

    errors = filter(f -> f.severity == :error, findings)
    warnings = filter(f -> f.severity == :warning, findings)

    println("review_config: $yaml_path")
    println("  $(length(errors)) error(s), $(length(warnings)) warning(s)")
    println("  (ported checks: dimensions/NPROC/CLOSEC, branch & grid geometry, JBDN, CST")
    println("  constituent activation -- ~1,150 more checks in pre_ivf_422.f90 not yet ported)")
    println()

    for f in errors
        println("ERROR   [$(f.section)] $(f.message)")
    end
    for f in warnings
        println("WARNING [$(f.section)] $(f.message)")
    end

    return length(errors)
end

if abspath(PROGRAM_FILE) == @__FILE__
    yaml_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..", "DetroitReservoir", "w2_config.yaml")
    n_errors = main(yaml_path)
    exit(n_errors > 0 ? 1 : 0)
end
