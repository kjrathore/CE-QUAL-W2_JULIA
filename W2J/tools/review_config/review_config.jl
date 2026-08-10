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
using Dates

struct Finding
    severity::Symbol   # :error or :warning
    section::String
    message::String
end

"Normalizes a config value that may be a bare scalar (NWB=1 case) or a YAML list (NWB>1)."
as_vector(x) = x isa AbstractVector ? x : Any[x]

"""
Upper-cased, stripped string form of a config value, or "" if the value is
YAML `null` (parsed by YAML.jl as Julia `nothing`) -- e.g. an unused PIPES
column when NPI=0. `string(nothing)` is the literal text "nothing", which
would otherwise silently pass through and fail an ON/OFF-style check as a
false positive. Every string-valued check below must go through this, not
a bare `string(x)`.
"""
str_or_empty(x) = x === nothing || x === missing ? "" : strip(uppercase(string(x)))

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

    closec = str_or_empty(dims["CLOSEC"])
    if !(closec in ("OFF", "ON"))
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
# Outlet structures -- pre_ivf_422.f90 lines 8896-8969
#
# GRID-DEPENDENT CHECKS SKIPPED: several source checks compare KTSTR/KBSTR/ESTR
# against EL(K,I) (layer elevation) and KB(I) (bottom active layer) -- both
# come from the geometry solve in init-geom.F90, not yet ported to W2J (see
# CLAUDE.md "Open questions"). Only the checks that need nothing but the
# control-file values themselves are ported here.
# ------------------------------------------------------------------------------
function check_structures(cfg)
    findings = Finding[]
    sec = find_section(cfg, "STRUCTURES for each branch")
    inout = find_section(cfg, "IN/OUTFLOW")
    grid = find_section(cfg, "BRANCH GRID")
    (sec === nothing || inout === nothing || grid === nothing) && return findings

    NST = inout["NST"]
    NSTR = as_vector(sec["NSTR"])
    DHS = as_vector(grid["DHS"])

    for jb in eachindex(NSTR)
        nstr_jb = NSTR[jb] isa Number ? NSTR[jb] : 0

        # pre_ivf_422.f90:8900-8903 -- DQ_EXTERNAL(JB) = DHS(JB)==0 (confirmed at
        # pre_ivf_422.f90:3084). No structures but an external downstream flow
        # boundary means there's nowhere for outflow to go.
        if jb <= length(DHS) && DHS[jb] == 0 && nstr_jb == 0
            push!(findings, Finding(:warning, "STRUCTURES", "branch $jb: no structures defined [NSTR=0] for external downstream flow boundary"))
        end
        # pre_ivf_422.f90:8965-8969 -- structures defined but no external flow boundary to use them
        if jb <= length(DHS) && DHS[jb] != 0 && nstr_jb > 0
            push!(findings, Finding(:error, "STRUCTURES", "branch $jb: NSTR=$nstr_jb > 0 but branch has an internal head/flow downstream boundary (DHS=$(DHS[jb]))"))
        end

        for js in 1:min(NST, 5)
            js > nstr_jb && break
            ktstr = as_vector(get(sec, "KTSTR$js", []))
            kbstr = as_vector(get(sec, "KBSTR$js", []))
            sinkc = as_vector(get(sec, "SINKC$js", []))
            wstr = as_vector(get(sec, "WSTR$js", []))
            jb > length(ktstr) && continue

            kt, kb = ktstr[jb], jb <= length(kbstr) ? kbstr[jb] : missing
            # pre_ivf_422.f90:8910-8912
            if kt isa Number && kt < 2
                push!(findings, Finding(:error, "STRUCTURES", "branch $jb structure $js: KTSTR=$kt < 2"))
            end
            # pre_ivf_422.f90:8914-8918
            if kt isa Number && kb isa Number && kt > kb
                push!(findings, Finding(:error, "STRUCTURES", "branch $jb structure $js: KTSTR=$kt > KBSTR=$kb"))
            end

            sk = jb <= length(sinkc) ? str_or_empty(sinkc[jb]) : ""
            # pre_ivf_422.f90:8939-8942
            if !(sk in ("LINE", "POINT"))
                push!(findings, Finding(:error, "STRUCTURES", "branch $jb structure $js: SINKC=\"$sk\" must be \"LINE\" or \"POINT\""))
            elseif sk == "LINE"
                # pre_ivf_422.f90:8944-8946
                w = jb <= length(wstr) ? wstr[jb] : missing
                if w isa Number && w <= 0.0
                    push!(findings, Finding(:error, "STRUCTURES", "branch $jb structure $js: WSTR=$w <= 0 while SINKC=LINE"))
                end
            end
        end
    end

    return findings
end

# ------------------------------------------------------------------------------
# Pipes -- pre_ivf_422.f90 lines 8976-9119
#
# GRID-DEPENDENT CHECKS SKIPPED: elevation-vs-EL(K,I)/KB(I) checks (EUPI,
# ETUPI/ETDPI, KBUPI/KBDPI vs KB) -- same reason as check_structures.
# ------------------------------------------------------------------------------
function check_pipes(cfg)
    findings = Finding[]
    sec = find_section(cfg, "PIPES")
    sec === nothing && return findings

    WPI = as_vector(get(sec, "WPI", []))
    DLXPI = as_vector(get(sec, "DLXPI", []))
    FPI = as_vector(get(sec, "FPI", []))
    FMINPI = as_vector(get(sec, "FMINPI", []))
    PUPIC = as_vector(get(sec, "PUPIC", []))
    PDPIC = as_vector(get(sec, "PDPIC", []))
    KTUPI = as_vector(get(sec, "KTUPI", []))
    KBUPI = as_vector(get(sec, "KBUPI", []))
    KTDPI = as_vector(get(sec, "KTDPI", []))
    KBDPI = as_vector(get(sec, "KBDPI", []))
    IDPI = as_vector(get(sec, "IDPI", []))
    LATPIC = as_vector(get(sec, "LATPIC", []))
    DYNPIPE = as_vector(get(sec, "DYNPIPE", []))

    n = maximum((length(WPI), length(DLXPI), length(FPI)); init=0)
    for jp in 1:n
        # pre_ivf_422.f90:9019-9026
        w = jp <= length(WPI) ? WPI[jp] : missing
        if w isa Number
            w <= 0.0 && push!(findings, Finding(:error, "PIPES", "pipe $jp: WPI=$w <= 0"))
            w > 5.0 && push!(findings, Finding(:warning, "PIPES", "pipe $jp: WPI=$w > 5 m"))
        end
        # pre_ivf_422.f90:9027-9034
        dlx = jp <= length(DLXPI) ? DLXPI[jp] : missing
        if dlx isa Number
            dlx <= 0.0 && push!(findings, Finding(:error, "PIPES", "pipe $jp: DLXPI=$dlx <= 0"))
            dlx > 1000.0 && push!(findings, Finding(:warning, "PIPES", "pipe $jp: DLXPI=$dlx > 1000 m"))
        end
        # pre_ivf_422.f90:9035-9038
        f = jp <= length(FPI) ? FPI[jp] : missing
        if f isa Number && (f <= 0.0 || f > 1.0)
            push!(findings, Finding(:error, "PIPES", "pipe $jp: FPI=$f must be in (0, 1]"))
        end
        # pre_ivf_422.f90:9039-9046
        fmin = jp <= length(FMINPI) ? FMINPI[jp] : missing
        if fmin isa Number
            fmin < 0.0 && push!(findings, Finding(:error, "PIPES", "pipe $jp: FMINPI=$fmin < 0"))
            fmin > 10.0 && push!(findings, Finding(:warning, "PIPES", "pipe $jp: FMINPI=$fmin > 10"))
        end
        # pre_ivf_422.f90:9047-9049
        pu = jp <= length(PUPIC) ? str_or_empty(PUPIC[jp]) : ""
        if !isempty(pu) && !(pu in ("SPECIFY", "DISTR", "DENSITY"))
            push!(findings, Finding(:error, "PIPES", "pipe $jp: PUPIC=\"$pu\" must be \"SPECIFY\", \"DISTR\", or \"DENSITY\""))
        end
        # pre_ivf_422.f90:9051-9064
        kt, kb = jp <= length(KTUPI) ? KTUPI[jp] : missing, jp <= length(KBUPI) ? KBUPI[jp] : missing
        if kt isa Number
            kt <= 1 && push!(findings, Finding(:error, "PIPES", "pipe $jp: KTUPI=$kt < 2"))
            if kb isa Number
                kt > kb && push!(findings, Finding(:error, "PIPES", "pipe $jp: KTUPI=$kt > KBUPI=$kb"))
                kt == kb && push!(findings, Finding(:warning, "PIPES", "pipe $jp: KTUPI=$kt == KBUPI=$kb"))
            end
        end
        # pre_ivf_422.f90:9065-9092,9105-9112 -- downstream/lateral fields only apply if IDPI != 0
        idpi = jp <= length(IDPI) ? IDPI[jp] : missing
        if idpi isa Number && idpi != 0
            pd = jp <= length(PDPIC) ? str_or_empty(PDPIC[jp]) : ""
            if !isempty(pd) && !(pd in ("SPECIFY", "DISTR", "DENSITY"))
                push!(findings, Finding(:error, "PIPES", "pipe $jp: PDPIC=\"$pd\" must be \"SPECIFY\", \"DISTR\", or \"DENSITY\""))
            end
            ktd, kbd = jp <= length(KTDPI) ? KTDPI[jp] : missing, jp <= length(KBDPI) ? KBDPI[jp] : missing
            if ktd isa Number
                ktd <= 1 && push!(findings, Finding(:error, "PIPES", "pipe $jp: KTDPI=$ktd < 2"))
                ktd isa Number && kbd isa Number && ktd > kbd && push!(findings, Finding(:error, "PIPES", "pipe $jp: KTDPI=$ktd > KBDPI=$kbd"))
            end
        end
        # pre_ivf_422.f90:9105-9108
        lat = jp <= length(LATPIC) ? str_or_empty(LATPIC[jp]) : ""
        if !isempty(lat) && !(lat in ("DOWN", "LAT"))
            push!(findings, Finding(:error, "PIPES", "pipe $jp: LATPIC=\"$lat\" must be \"DOWN\" or \"LAT\""))
        end
        # pre_ivf_422.f90:9109-9112
        dyn = jp <= length(DYNPIPE) ? str_or_empty(DYNPIPE[jp]) : ""
        if !isempty(dyn) && !(dyn in ("ON", "OFF"))
            push!(findings, Finding(:error, "PIPES", "pipe $jp: DYNPIPE=\"$dyn\" must be \"ON\" or \"OFF\""))
        end
    end

    return findings
end

# ------------------------------------------------------------------------------
# Withdrawals -- pre_ivf_422.f90 lines 9931-9968
#
# GRID-DEPENDENT CHECKS SKIPPED: EWD vs EL(2,I)/EL(KB+1,I), KBWD vs KB(IWD) --
# same reason as check_structures.
# ------------------------------------------------------------------------------
function check_withdrawals(cfg)
    findings = Finding[]
    sec = find_section(cfg, "WITHDRAWALS")
    grid = find_section(cfg, "BRANCH GRID")
    (sec === nothing || grid === nothing) && return findings

    IWD = as_vector(get(sec, "IWD", []))
    KTWD = as_vector(get(sec, "KTWD", []))
    KBWD = as_vector(get(sec, "KBWD", []))
    US = as_vector(grid["US"]); DS = as_vector(grid["DS"])

    for jw in eachindex(IWD)
        iwd = IWD[jw]
        iwd isa Number || continue

        # pre_ivf_422.f90:9940-9947
        for jb in eachindex(US)
            if iwd == US[jb] - 1 || (jb == 1 && iwd == 0)
                push!(findings, Finding(:error, "WITHDRAWALS", "withdrawal $jw: IWD=$iwd is an upstream boundary segment in branch $jb"))
            elseif iwd == DS[jb] + 1
                push!(findings, Finding(:error, "WITHDRAWALS", "withdrawal $jw: IWD=$iwd is a downstream boundary segment in branch $jb"))
            end
        end

        # pre_ivf_422.f90:9949-9968
        kt = jw <= length(KTWD) ? KTWD[jw] : missing
        kb = jw <= length(KBWD) ? KBWD[jw] : missing
        if kt isa Number
            kt < 2 && push!(findings, Finding(:error, "WITHDRAWALS", "withdrawal $jw: KTWD=$kt < 2"))
            kt isa Number && kb isa Number && kt > kb && push!(findings, Finding(:error, "WITHDRAWALS", "withdrawal $jw: KTWD=$kt > KBWD=$kb"))
        end
    end

    return findings
end

# ------------------------------------------------------------------------------
# Tributaries -- pre_ivf_422.f90 lines 9976-10010
#
# GRID-DEPENDENT CHECKS SKIPPED: ETTR/EBTR vs EL(2,I)/EL(KB+1,I) -- same
# reason as check_structures. The ETTR<EBTR ordering check needs no grid data
# so it IS ported.
# ------------------------------------------------------------------------------
function check_tributaries(cfg)
    findings = Finding[]
    sec = find_section(cfg, "TRIB PLACEMENT")
    grid = find_section(cfg, "BRANCH GRID")
    (sec === nothing || grid === nothing) && return findings

    PTRC = as_vector(get(sec, "PTRC", []))
    ITR = as_vector(get(sec, "ITR", []))
    ELTRT = as_vector(get(sec, "ELTRT", []))
    ELTRB = as_vector(get(sec, "ELTRB", []))
    US = as_vector(grid["US"]); DS = as_vector(grid["DS"])

    for jt in eachindex(PTRC)
        # pre_ivf_422.f90:9976-9979
        ptrc = str_or_empty(PTRC[jt])
        isempty(ptrc) && continue
        if !(ptrc in ("SPECIFY", "DISTR", "DENSITY"))
            push!(findings, Finding(:error, "TRIB PLACEMENT", "tributary $jt: PTRC=\"$ptrc\" must be \"SPECIFY\", \"DISTR\", or \"DENSITY\""))
        end

        # pre_ivf_422.f90:9980-9987
        itr = jt <= length(ITR) ? ITR[jt] : missing
        if itr isa Number && !any(jb -> itr >= US[jb] && itr <= DS[jb], eachindex(US))
            push!(findings, Finding(:error, "TRIB PLACEMENT", "tributary $jt: ITR=$itr is a boundary segment or not in the active grid"))
        end

        # pre_ivf_422.f90:9988,10005-10009 -- elevation ordering only checked when PTRC=SPECIFY
        if ptrc == "SPECIFY"
            top = jt <= length(ELTRT) ? ELTRT[jt] : missing
            bot = jt <= length(ELTRB) ? ELTRB[jt] : missing
            if top isa Number && bot isa Number && top < bot
                push!(findings, Finding(:error, "TRIB PLACEMENT", "tributary $jt: ELTRT=$top < ELTRB=$bot"))
            end
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
        active = str_or_empty(get(c, "active", ""))
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
            pv = str_or_empty(v)
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
    append!(findings, check_structures(cfg))
    append!(findings, check_pipes(cfg))
    append!(findings, check_withdrawals(cfg))
    append!(findings, check_tributaries(cfg))
    append!(findings, check_cst(cfg))
    return findings
end

"Builds the full report as a single string, so it prints and writes to `Outputs/reviewed.log` identically."
function report_text(yaml_path, findings)
    errors = filter(f -> f.severity == :error, findings)
    warnings = filter(f -> f.severity == :warning, findings)

    io = IOBuffer()
    println(io, "review_config: $yaml_path")
    println(io, "run at: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "  $(length(errors)) error(s), $(length(warnings)) warning(s)")
    println(io, "  (ported checks: dimensions/NPROC/CLOSEC, branch & grid geometry, JBDN, structures,")
    println(io, "  pipes, withdrawals, tributaries, CST constituent activation -- grid-elevation-")
    println(io, "  dependent checks (need EL(K,I)/KB(I) from init-geom.F90, not yet ported to W2J)")
    println(io, "  are skipped; ~1,100 more checks in pre_ivf_422.f90 not yet ported at all)")
    println(io)

    for f in errors
        println(io, "ERROR   [$(f.section)] $(f.message)")
    end
    for f in warnings
        println(io, "WARNING [$(f.section)] $(f.message)")
    end

    return String(take!(io))
end

"""
    main(yaml_path; log_path=nothing)

Runs the checks and prints the report to stdout, same as before. Also writes
the identical report to a log file so it's a durable, revisitable artifact --
matching the project convention of writing generated output under
`DetroitReservoir/Outputs/` (see `Plotting/LongitudinalProfile.jl`'s test).
Defaults to `<yaml_path's directory>/Outputs/reviewed.log`; pass `log_path`
to write somewhere else.
"""
function main(yaml_path::AbstractString; log_path::Union{AbstractString,Nothing}=nothing)
    cfg = YAML.load_file(yaml_path)
    findings = run_checks(cfg)
    text = report_text(yaml_path, findings)

    print(text)

    log_path = log_path === nothing ? joinpath(dirname(yaml_path), "Outputs", "reviewed.log") : log_path
    mkpath(dirname(log_path))
    open(log_path, "w") do f
        write(f, text)
    end
    println("(report written to $(abspath(log_path)))")

    return length(filter(f -> f.severity == :error, findings))
end

if abspath(PROGRAM_FILE) == @__FILE__
    yaml_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..", "DetroitReservoir", "w2_config.yaml")
    log_path = length(ARGS) >= 2 ? ARGS[2] : nothing
    n_errors = main(yaml_path; log_path=log_path)
    exit(n_errors > 0 ? 1 : 0)
end
