# ==============================================================================
# tools/xlsx_to_yaml/convert.jl
#
# Converts a CE-QUAL-W2 "w2_con" Excel workbook (the .xlsm master, which is the
# real source of truth -- w2_con.csv is exported FROM it) into a single
# w2_config.yaml. This is a standalone conversion utility, not part of the W2J
# package itself -- kept in its own tools/ subfolder with its own Project.toml
# (just XLSX.jl) so the runtime package doesn't carry a spreadsheet-reading
# dependency.
#
# WHY YAML: w2_con.csv's row-position-only format (no self-describing labels
# reaching InputReader.jl -- see that file's Phase-A/Tier-1 scope notes) forces
# every future parser change to re-derive exact row counts by hand. The Excel
# workbook, by contrast, already carries real structure -- column B holds a
# section title or a "VARNAME - description" label, column C onward holds
# either data or column headers. This script exploits that structure to emit a
# labeled, commented, human-editable YAML file. Future InputReader work can
# target this YAML instead of re-deriving w2_con.csv's row offsets.
#
# SHEET LAYOUT (confirmed against the Detroit workbook's "w2_con.csv" sheet,
# 942 rows x 56 cols):
#   - Column A: free-standing commentary, NOT aligned with the same row's
#     data -- explicitly flagged by the sheet's own row 1 ("Note COL A and B
#     are not written out to w2_con.csv"). Ignored entirely by this script.
#   - Column B: either a section title on a header row, or "VARNAME -
#     description" / "VARNAME description" on a field row.
#   - Column C onward: on a header row, either (a) generic per-instance
#     placeholders ("WB1","WB2",...,"BR1",...,"PIPE1",...) -- an entity
#     block, where every following field row holds one value per instance
#     column, (b) a repeated single token ("DLTD","DLTD",...) -- a Fortran
#     array, whose value(s) sit in the next row, or (c) anything else --
#     resolved by lookahead (see below), not by guessing at the header shape.
#
# CONFIRMED INVARIANT (checked against every section header in the Detroit
# workbook): a header row is ALWAYS immediately preceded by a blank-separator
# row (column B missing, column C a lone space " "), and a block's field rows
# always run until the NEXT blank separator. This is the only thing used to
# decide where a block ends -- earlier drafts of this script tried to
# recognize "does this row look like a new header" from row shape alone, and
# that broke on ordinary field rows whose *value* happened to look like a
# header token (e.g. a field literally valued "ON" looks identical, in
# isolation, to a one-token dimension-header). Don't reintroduce that check.
#
# Once a header's block is scanned, two shapes remain for the general (non-
# entity, non-array) case, disambiguated by how many rows sit in the block
# before the next blank separator:
#   - Exactly one row: a dimension-scalar block -- the header names N fields,
#     that one row supplies N values 1:1 (e.g. "NWB, NBR, IMX, KMX, NPROC,
#     CLOSEC" -> "1, 4, 31, 117, 1, OFF").
#   - More than one row: a vertical field block -- each row supplies its OWN
#     key from column B (split on the first " - "/"-"/"#" separator) plus its
#     own value(s) from column C onward. This also correctly handles headers
#     whose token list isn't a WBn/BRn placeholder run but is still just
#     labeling value columns (e.g. "HYD PRINT"'s "HNAME, FMTH" header, whose
#     ~16 following rows are U, W, T, RHO, ... each supplying its own name).
#
# KNOWN LIMITATION: the constituent tables (CST / CST DERI, ~65 rows total)
# are wide tables whose header names real per-column attributes (CNAME2,
# CNAME, CAC, FMTC, ...) and whose data rows carry a short constituent code
# in column B (e.g. "TDS") rather than a "VARNAME - description" label. This
# script's vertical-field fallback still applies (key = the constituent code,
# values = the rest of that row), which is inspectable and round-trippable
# but loses the header's own column names as explicit per-value keys. Fixing
# this needs a dedicated special case (match on "CST" in the header text,
# shift the value list left by one to fold column B into the aligned data) --
# flagged here rather than guessed at, per this project's working discipline.
#
# NOT PERFECT BY DESIGN otherwise, too: ~900 rows of a 25-year-old Fortran
# control file have more format quirks than any fixed ruleset catches on the
# first pass. Anything this script can't confidently classify is preserved
# verbatim under `_raw_unclassified` (row number + full row contents) instead
# of being dropped or mis-mapped. Spot-check new sections against the source
# workbook before trusting them.
# ==============================================================================

using XLSX

const SHEET_NAME = "w2_con.csv"

# ------------------------------------------------------------------------------
# Cell / row helpers
# ------------------------------------------------------------------------------

"Trim trailing `missing` entries off a row slice."
function trimmed(row::AbstractVector)
    last_real = findlast(x -> x !== missing, row)
    return last_real === nothing ? eltype(row)[] : row[1:last_real]
end

is_blank_marker(b, rest) = b === missing && (isempty(rest) || (rest[1] isa AbstractString && strip(rest[1]) == ""))

"""
A block ends at a real blank separator OR at a row whose column B is neither
missing nor a string -- e.g. row 400 in the Detroit workbook, a stray bare
integer (`40`) left over from an Excel formula/annotation ("Your last column
should be in COLUMN: N"), with no string label and no useful data alongside
it. Every genuine field row in every confirmed block has a string column B;
treating a bare-number column B as "still part of the block" was silently
merging the CST constituent table into whatever entity block preceded it.
"""
ends_block(b, rest) = is_blank_marker(b, rest) || (b !== missing && !(b isa AbstractString))

"Generic per-instance token, e.g. \"WB1\", \" BR2\", \"PIPE10\", \"MacGroup3\"."
function entity_match(s)
    s isa AbstractString || return nothing
    m = match(r"^\s*([A-Za-z_]+)\s*([0-9]{1,3})\s*$", s)
    return m
end

"Is `tokens` a run of generic entity placeholders (WB1, WB2, ... or BR1, BR2, ...)?"
function looks_like_entity_header(tokens)
    length(tokens) < 1 && return false
    ms = [entity_match(t) for t in tokens]
    any(isnothing, ms) && return false
    return length(unique(m -> uppercase(strip(m.captures[1])), ms)) == 1
end

"Is `tokens` all the same literal string repeated (a Fortran array header, e.g. DLTD DLTD DLTD)?"
function looks_like_repeated_array(tokens)
    length(tokens) < 2 && return false
    all(t -> t isa AbstractString, tokens) || return false
    normed = uppercase.(strip.(tokens))
    return length(unique(normed)) == 1
end

"Extract (key, comment) from a field-row label like \"VBC - volume balance computation\" or \"NLMIN # of layers\"."
function split_label(label::AbstractString)
    s = strip(label)
    m = match(r"^([A-Za-z][A-Za-z0-9_%]*)\s*[-#,:]?\s*(.*)$", s)
    if m === nothing
        return (s, "")
    end
    return (uppercase(strip(m.captures[1])), strip(m.captures[2]))
end

# ------------------------------------------------------------------------------
# YAML text emission (hand-rolled, not via a generic serializer, so we keep
# full control over inline `#` comments -- YAML.jl's writer does not support
# round-tripping comments, and comments are the whole point of this format
# per the project decision to move off w2_con.csv's opaque row-positional
# format).
# ------------------------------------------------------------------------------

using Dates

"""
Single-quoted YAML scalar: only `'` needs escaping (doubled), so this is safe
for Windows paths (backslashes are NOT an escape character in single-quoted
YAML, unlike double-quoted) and for values that already carry embedded
double quotes (several Excel cells do, e.g. Fortran-format strings).
"""
function yaml_scalar(v)
    v === missing && return "null"
    if v isa AbstractString
        s = replace(strip(v), "'" => "''")
        return "'$s'"
    elseif v isa Number
        return string(v)
    elseif v isa Dates.DateTime || v isa Dates.Date
        return "'$(string(v))'"
    else
        return "'$(string(v))'"
    end
end

yaml_key(k::AbstractString) = occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", k) ? k : "\"$(replace(k, "\"" => "\\\""))\""

emit(io, indent, s) = println(io, " "^indent * s)

function emit_field!(io, indent, key, values, comment)
    ck = yaml_key(key)
    line = if length(values) == 1
        "$ck: $(yaml_scalar(values[1]))"
    else
        "$ck: [" * join(yaml_scalar.(values), ", ") * "]"
    end
    if !isempty(comment)
        line *= "  # $(comment)"
    end
    emit(io, indent, line)
end

"""
    emit_cst_block!(io, indent, block_rows, data, is_deri, nwb, nbr, ntr)

Special case for the "CST - Concentration State variables and initial
conditions" and "CST DERI - Derived concentration state variables" tables.
These don't fit the generic vertical-field shape (column B on their data
rows holds a bare row index for CST, or the constituent's own short name for
CST DERI -- neither is a "VARNAME - description" label), and their column
count depends on NWB/NBR/NTR (captured earlier in the same pass, from the
"GRID/NPROC/CLOSE DIALOG BOX" and "IN/OUTFLOW" dimension blocks -- see the
manual, `W2manual455_Part3_InputOutputFiles_rev1.pdf` pp.156-157 and p.139,
for the confirmed column order this mirrors).

CST column order after CNAME2/CNAME/CAC/FMTC/CMULT: C2IWB, CPRWBC, and
C_Atm_Deposition each repeat once per waterbody (NWB); CINBRC repeats once
per branch (NBR); CTRTRC repeats once per tributary (NTR, minimum 1); CDTBRC
and CPRBRC each repeat once per branch (NBR) again.

CST DERI column order after CDNAME2/CDNAME/FMTCD/CDMULT: CDWBC repeats once
per waterbody (NWB).
"""
function emit_cst_block!(io, indent, block_rows, data, is_deri, nwb, nbr, ntr)
    ind = " "^indent
    for row in block_rows
        fb = data[row, 2]
        frest = trimmed(data[row, 3:end])
        isempty(frest) && continue

        if is_deri
            # Column B is blank here -- CDNAME2 (the short name) is the first
            # column-C token, same alignment as the header row (unlike plain
            # CST, where column B holds a bare row index instead).
            shortname = frest[1] isa AbstractString ? strip(frest[1]) : string(frest[1])
            isempty(shortname) && continue
            longname = length(frest) >= 2 ? frest[2] : missing
            fmt = length(frest) >= 3 ? frest[3] : missing
            mult = length(frest) >= 4 ? frest[4] : missing
            printwb = length(frest) >= 5 ? frest[5:end] : Any[]
            nwb > 0 && (printwb = printwb[1:min(nwb, length(printwb))])
            println(io, "$(ind)$(yaml_key(shortname)):")
            println(io, "$(ind)  long_name: $(yaml_scalar(longname))")
            println(io, "$(ind)  format: $(yaml_scalar(fmt))")
            println(io, "$(ind)  multiplier: $(yaml_scalar(mult))")
            println(io, "$(ind)  print_wb: [" * join(yaml_scalar.(printwb), ", ") * "]")
        else
            # Column B is a bare row index here, not a field identity -- ignored.
            shortname = frest[1] isa AbstractString ? strip(frest[1]) : string(frest[1])
            isempty(shortname) && continue
            idx = Ref(2)
            take1() = (v = idx[] <= length(frest) ? frest[idx[]] : missing; idx[] += 1; v)
            taken(n) = begin
                lo, hi = idx[], min(idx[] + n - 1, length(frest))
                vals = hi >= lo ? frest[lo:hi] : Any[]
                idx[] += n
                vals
            end
            longname = take1()
            active = take1()
            fmt = take1()
            mult = take1()
            c2iwb = taken(nwb)
            cprwbc = taken(nwb)
            atm = taken(nwb)
            cinbrc = taken(nbr)
            ctrtrc = taken(max(ntr, 1))
            cdtbrc = taken(nbr)
            cprbrc = taken(nbr)
            println(io, "$(ind)$(yaml_key(shortname)):")
            println(io, "$(ind)  long_name: $(yaml_scalar(longname))")
            println(io, "$(ind)  active: $(yaml_scalar(active))")
            println(io, "$(ind)  format: $(yaml_scalar(fmt))")
            println(io, "$(ind)  multiplier: $(yaml_scalar(mult))")
            println(io, "$(ind)  initial_conc: [" * join(yaml_scalar.(c2iwb), ", ") * "]  # per waterbody (NWB=$nwb)")
            println(io, "$(ind)  print_snp: [" * join(yaml_scalar.(cprwbc), ", ") * "]  # per waterbody")
            println(io, "$(ind)  atm_deposition: [" * join(yaml_scalar.(atm), ", ") * "]  # per waterbody")
            println(io, "$(ind)  inflow: [" * join(yaml_scalar.(cinbrc), ", ") * "]  # per branch (NBR=$nbr)")
            println(io, "$(ind)  tributary: [" * join(yaml_scalar.(ctrtrc), ", ") * "]  # per tributary (NTR=$ntr)")
            println(io, "$(ind)  distributed: [" * join(yaml_scalar.(cdtbrc), ", ") * "]  # per branch")
            println(io, "$(ind)  precipitation: [" * join(yaml_scalar.(cprbrc), ", ") * "]  # per branch")
        end
    end
end

"""
    emit_vertical_fields!(io, indent, rows, data, raw_unclassified)

Emits one field per row (key from column B via `split_label`, value(s) from
column C onward) -- EXCEPT when consecutive rows share the same extracted
key (e.g. 5 macrophyte-group placeholder rows all labeled "MAC Waterbody
macrophyte <n> computations..." -- the workbook describes the group number
in free text rather than a MAC1/MAC2-style suffix, so `split_label` alone
can't tell them apart). Those runs are merged into one array field instead
of being silently overwritten under a duplicate YAML key.
"""
function emit_vertical_fields!(io, indent, rows, data, raw_unclassified)
    parsed = Tuple{String,String,Vector}[]
    for row in rows
        fb = data[row, 2]
        frest = trimmed(data[row, 3:end])
        if fb isa AbstractString
            key, comment = split_label(fb)
            push!(parsed, (key, comment, isempty(frest) ? Any[missing] : frest))
        else
            push!(raw_unclassified, (row, trimmed(data[row, 2:end])))
        end
    end

    i = 1
    while i <= length(parsed)
        key, comment, vals = parsed[i]
        j = i + 1
        while j <= length(parsed) && parsed[j][1] == key
            vals = vcat(vals, parsed[j][3])
            j += 1
        end
        emit_field!(io, indent, key, vals, comment)
        i = j
    end
end

# ------------------------------------------------------------------------------
# Main conversion
# ------------------------------------------------------------------------------

function convert_xlsx_to_yaml(xlsx_path::AbstractString, yaml_path::AbstractString)
    xf = XLSX.readxlsx(xlsx_path)
    sh = xf[SHEET_NAME]
    data = sh[:]
    nrows, ncols = size(data)

    io = IOBuffer()
    println(io, "# Auto-generated by W2J/tools/xlsx_to_yaml/convert.jl")
    println(io, "# Source workbook: $(basename(xlsx_path)), sheet \"$SHEET_NAME\"")
    println(io, "# Best-effort structural transcription -- see the module docstring in")
    println(io, "# convert.jl for the row-classification rules and their limits. Rows the")
    println(io, "# converter could not confidently classify are under `_raw_unclassified`.")
    println(io)

    raw_unclassified = Vector{Tuple{Int,Vector}}()

    # Populated from the "GRID/NPROC/CLOSE DIALOG BOX" and "IN/OUTFLOW" dimension
    # blocks (both precede the constituent tables in every known layout of this
    # file) -- needed to correctly slice the CST / CST DERI tables' per-waterbody /
    # per-branch / per-tributary columns. See emit_cst_block!'s docstring.
    nwb = Ref(0)
    nbr = Ref(0)
    ntr = Ref(0)

    # Skip the fixed version-banner preamble (rows 1-2: "w2_con.csv file format" /
    # "Control File version", CE-QUAL-W2 version number) -- workbook boilerplate,
    # not model configuration, and not shaped like any of the row patterns below.
    # Start real processing at the "TITLE C" marker row.
    r = 1
    while r <= nrows && !(data[r, 2] isa AbstractString && strip(data[r, 2]) == "TITLE C")
        r += 1
    end

    while r <= nrows
        b = data[r, 2]
        rest = trimmed(data[r, 3:end])

        if is_blank_marker(b, rest)
            r += 1
            continue
        end

        # --- Title block: "TITLE C" marker -> next 10 rows' column C are raw text lines ---
        if b isa AbstractString && strip(b) == "TITLE C"
            println(io, "title:")
            for tr in (r+1):(r+10)
                tr > nrows && break
                tval = data[tr, 3]
                tval === missing && continue
                emit(io, 2, "- $(yaml_scalar(tval))")
            end
            println(io)
            r += 11
            continue
        end

        # --- Header row candidates (only ever reached right after a blank separator) ---
        if b isa AbstractString && !isempty(rest)
            section_title = strip(b)

            if looks_like_repeated_array(rest)
                # Fortran array: name repeated across header cells, value(s) on next row.
                arr_name = uppercase(strip(rest[1]))
                vr = r + 1
                vals = vr <= nrows ? trimmed(data[vr, 3:end]) : []
                println(io, "$(yaml_key(section_title)):")
                emit_field!(io, 2, arr_name, isempty(vals) ? [missing] : vals, "")
                println(io)
                r = vr + 1
                continue

            elseif looks_like_entity_header(rest)
                # Entity block (WB1/BR1/PIPE1/...): every row until the next blank
                # separator is one field, keyed off its own column B.
                fr = r + 1
                while fr <= nrows
                    fb = data[fr, 2]
                    frest = trimmed(data[fr, 3:end])
                    ends_block(fb, frest) && break
                    fr += 1
                end
                println(io, "$(yaml_key(section_title)):")
                emit_vertical_fields!(io, 2, (r+1):(fr-1), data, raw_unclassified)
                println(io)
                r = fr
                continue

            else
                # General case: look ahead to the next blank separator to decide
                # dimension-scalar (exactly one data row) vs. vertical-field
                # (more than one) -- see module docstring for why this is decided
                # by lookahead rather than by the header row's own shape.
                # NOTE: uses is_blank_marker, not ends_block -- CST/CST DERI data
                # rows legitimately have a non-string column B (a bare row index
                # for CST, see emit_cst_block!), so this loop can't treat that as
                # a block terminator the way the entity-block loop above does.
                block_rows = Int[]
                tr = r + 1
                while tr <= nrows
                    tb = data[tr, 2]
                    trest = trimmed(data[tr, 3:end])
                    is_blank_marker(tb, trest) && break
                    push!(block_rows, tr)
                    tr += 1
                end

                println(io, "$(yaml_key(section_title)):")
                if occursin("CST", uppercase(section_title))
                    emit_cst_block!(io, 2, block_rows, data, occursin("DERI", uppercase(section_title)), nwb[], nbr[], ntr[])
                elseif length(block_rows) <= 1
                    vals = isempty(block_rows) ? Any[] : trimmed(data[block_rows[1], 3:end])
                    for (i, tok) in enumerate(rest)
                        key = uppercase(strip(tok isa AbstractString ? tok : string(tok)))
                        v = i <= length(vals) ? vals[i] : missing
                        emit_field!(io, 2, key, [v], "")
                        key == "NWB" && v isa Number && (nwb[] = Int(v))
                        key == "NBR" && v isa Number && (nbr[] = Int(v))
                        key == "NTR" && v isa Number && (ntr[] = Int(v))
                    end
                else
                    emit_vertical_fields!(io, 2, block_rows, data, raw_unclassified)
                end
                println(io)
                r = isempty(block_rows) ? r + 1 : block_rows[end] + 1
                continue
            end
        end

        # --- Anything else: flag for manual review, don't guess. ---
        push!(raw_unclassified, (r, trimmed(data[r, 2:end])))
        r += 1
    end

    if !isempty(raw_unclassified)
        println(io, "_raw_unclassified:")
        for (row, vals) in raw_unclassified
            valstr = join(yaml_scalar.(vals), ", ")
            emit(io, 2, "- {row: $row, cells: [$valstr]}")
        end
    end

    open(yaml_path, "w") do f
        write(f, String(take!(io)))
    end

    return (n_raw = length(raw_unclassified), n_rows = nrows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    xlsx_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..", "DetroitReservoir", "w2_con_DetroitReservoir4.5.xlsm")
    yaml_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "..", "..", "DetroitReservoir", "w2_config.yaml")
    result = convert_xlsx_to_yaml(xlsx_path, yaml_path)
    println("Converted $(xlsx_path) -> $(yaml_path)")
    println("$(result.n_rows) rows scanned, $(result.n_raw) rows unclassified (see _raw_unclassified in the output).")
end
