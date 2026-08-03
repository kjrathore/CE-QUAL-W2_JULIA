# ==============================================================================
# Core/State.jl
#
# Foundation shared-state structs for the W2J port.
# Source: w2modules.F90, MODULE GLOBAL / GEOMC / NAMESC (the universal `USE`
# dependency in the original Fortran codebase -- everyone depends on this,
# it depends on nothing, and it calls nothing itself).
#
# DESIGN RULE (per project decision log #3 / #2): field names are kept
# IDENTICAL to the Fortran source (NWB, NBR, IMX, KMX, ...) so the existing
# W2 user manual, control-file labels, and decades of literature cross-
# reference directly with zero translation. The w2_con.csv label lines
# already use these exact tokens -- confirmed against w2modules.F90.
#
# NOT YET PORTED HERE (deliberately):
#   - The actual w2_con.csv / input.F90 read logic. w2modules.F90 only
#     DECLARES storage; the real field read ORDER lives in input.F90
#     (Tier 2), which has not been provided yet. Building the parser
#     before seeing that file risks a silent field-order mismatch.
#   - STRUCTURES, TRANS, SURFHE, TVDC, KINETIC module data -- these belong
#     to Hydrodynamics/Structures.jl, Hydrodynamics/Transport.jl, etc. per
#     the package map in W2J_README.md, not to the foundation struct.
#
# OPEN DESIGN NOTE: the Fortran GLOBAL module declares U, W, T2, AZ, RHO,
# ST, SB, DLTLIM, VSH, ADMX, DM, ADMZ, HDG, HPG, GRAV as POINTER (not
# ALLOCATABLE). This is almost certainly the old-step/new-step array-swap
# optimization (avoid copying T1<->T2 each timestep by swapping pointers
# instead). Julia needs no equivalent trick -- swapping which array is
# "current" is just rebinding a name, free either way -- so these are
# deliberately NOT struct fields yet. Confirm against temperature.F90 /
# the actual swap call sites before deciding the final representation
# (likely: just two plain arrays, swapped by the timestep loop itself).
# ==============================================================================

module W2Core

export W2Global, W2Geometry, W2Names, W2TimeControl

# ------------------------------------------------------------------------------
# MODULE GLOBAL
# ------------------------------------------------------------------------------
"""
    W2Global

Mirrors `MODULE GLOBAL` in w2modules.F90.

Grid-defining scalars (NWB, NBR, IMX, KMX, NPROC, CLOSEC) correspond exactly
to the w2_con.csv header line `NWB, NBR, IMX, KMX, NPROC, CLOSEC` -- confirmed
field-for-field against the uploaded control file.

All other ALLOCATABLE arrays start empty (size-0) and are sized once
NWB/NBR/IMX/KMX are known, mirroring the Fortran two-phase
read-dimensions-then-allocate pattern (the allocation step itself belongs
to IO/InputReader.jl, not here).
"""
mutable struct W2Global
    # --- grid-defining scalars : w2_con.csv "NWB, NBR, IMX, KMX, NPROC, CLOSEC" ---
    NWB::Int
    NBR::Int
    IMX::Int
    KMX::Int
    NPROC::Int      # Fortran comment: "# of processors (INACTIVE at this time)".
                    # REPURPOSED in W2J: read from the same w2_con.csv field and
                    # used to bound the thread count for Threads.@threads loops
                    # over segments/cells. No new config surface needed.
    CLOSEC::Bool    # ON/OFF string in CSV -> Bool

    # --- w2_con.csv "NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU" ---
    NTR::Int; NST::Int; NIW::Int; NWD::Int
    NGT::Int; NSP::Int; NPI::Int; NPU::Int

    # --- w2_con.csv "NGC, NSS, NAL, NEP, NBOD, NMC, NZP" ---
    NGC::Int; NSS::Int; NAL::Int; NEP::Int
    NBOD::Int; NMC::Int; NZP::Int

    # --- derived / bookkeeping counts (not in w2_con.csv directly) ---
    NCT::Int; NTR1::Int
    NWDO::Int; NIKTSR::Int; NUNIT::Int
    NOD::Int
    NDC::Int          # Fortran default: NDC = 27
    NHY::Int          # Fortran default: NHY = 15
    NFL::Int          # Fortran default: NFL = 142
    NEPT::Int; NZPT::Int; NMCT::Int
    NGCS::Int; NGCE::Int
    NZOOS::Int; NZOOE::Int

    # --- timestep control ---
    DLT::Float64; DLTMIN::Float64; DLTTVD::Float64; DZMAX::Float64
    BETABR::Float64; START::Float64; HMAX2::Float64; CURRENT::Float64

    # --- per-segment / per-layer working arrays (TARGET, ALLOCATABLE in Fortran) ---
    T1::Matrix{Float64}; TSS::Matrix{Float64}
    C1::Array{Float64,3}; C2::Array{Float64,3}; C1S::Array{Float64,3}
    CSSB::Array{Float64,3}; CSSK::Array{Float64,3}
    KF::Array{Float64,3}; CD::Array{Float64,3}
    HYD::Array{Float64,3}
    AF::Array{Float64,4}; EF::Array{Float64,4}

    ICETH::Vector{Float64}; ELKT::Vector{Float64}; HMULT::Vector{Float64}
    CMULT::Vector{Float64}; CDMULT::Vector{Float64}; WIND2::Vector{Float64}
    AZMAX::Vector{Float64}; PALT::Vector{Float64}; Z0::Vector{Float64}

    QSS::Matrix{Float64}; VOLUH2::Matrix{Float64}; VOLDH2::Matrix{Float64}
    QUH1::Matrix{Float64}; QDH1::Matrix{Float64}
    UXBR::Matrix{Float64}; UYBR::Matrix{Float64}; VOL::Matrix{Float64}

    ALLIM::Array{Float64,3}; APLIM::Array{Float64,3}
    ANLIM::Array{Float64,3}; ASLIM::Array{Float64,3}; KFS::Array{Float64,3}
    ELLIM::Array{Float64,3}; EPLIM::Array{Float64,3}
    ENLIM::Array{Float64,3}; ESLIM::Array{Float64,3}

    # --- connectivity / indexing (INTEGER ALLOCATABLE in Fortran) ---
    BS::Vector{Int}; BE::Vector{Int}; US::Vector{Int}; CUS::Vector{Int}
    DS::Vector{Int}; JBDN::Vector{Int}
    KB::Vector{Int}; KTI::Vector{Int}; SKTI::Vector{Int}; KTWB::Vector{Int}
    KBMIN::Vector{Int}; CDHS::Vector{Int}
    UHS::Vector{Int}; DHS::Vector{Int}; UQB::Vector{Int}; DQB::Vector{Int}
    OPT::Matrix{Int}
    NBODC::Vector{Int}; NBODN::Vector{Int}; NBODP::Vector{Int}

    ICE::Vector{Bool}; ICE_CALC::Vector{Bool}; LAYERCHANGE::Vector{Bool}
    BR_INACTIVE::Vector{Bool}; BR_NOTECPLOT::Vector{Bool}

    # --- misc ---
    RSIFN::String; MODDIR::String
end

"""
    W2Global()

Default-constructs an empty `W2Global` with zeroed scalars and zero-length
arrays. Real sizing happens once NWB/NBR/IMX/KMX are read from w2_con.csv
(see IO/InputReader.jl, not yet implemented -- pending input.F90).
"""
function W2Global()
    W2Global(
        0, 0, 0, 0, 1, false,                     # NWB..CLOSEC
        0, 0, 0, 0, 0, 0, 0, 0,                   # NTR..NPU
        0, 0, 0, 0, 0, 0, 0,                      # NGC..NZP
        0, 0, 0, 0, 0,                            # NCT, NTR1, NWDO, NIKTSR, NUNIT
        0, 27, 15, 142,                           # NOD, NDC, NHY, NFL
        0, 0, 0, 0, 0, 0, 0,                      # NEPT..NZOOE
        0.0, 0.0, 0.0, 100.0,                     # DLT..DZMAX (Fortran default DZMAX=1.0D2)
        0.0, 0.0, 0.0, 0.0,                       # BETABR..CURRENT
        zeros(0, 0), zeros(0, 0),                 # T1, TSS
        zeros(0, 0, 0), zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0),
        zeros(0, 0, 0, 0), zeros(0, 0, 0, 0),
        Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[],
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0),
        zeros(0, 0, 0), zeros(0, 0, 0),
        Int[], Int[], Int[], Int[],
        Int[], Int[],
        Int[], Int[], Int[], Int[],
        Int[], Int[],
        Int[], Int[], Int[], Int[],
        zeros(Int, 0, 0),
        Int[], Int[], Int[],
        Bool[], Bool[], Bool[],
        Bool[], Bool[],
        "", "",
    )
end

# ------------------------------------------------------------------------------
# MODULE GEOMC
# ------------------------------------------------------------------------------
"""
    W2Geometry

Mirrors `MODULE GEOMC` -- per-segment / per-layer geometry (widths, depths,
slopes, bed elevations). Dimensioned by IMX/KMX once read.
"""
mutable struct W2Geometry
    JBUH::Vector{Int}; JBDH::Vector{Int}; JWUH::Vector{Int}; JWDH::Vector{Int}

    ALPHA::Vector{Float64}; SINA::Vector{Float64}; COSA::Vector{Float64}
    SLOPE::Vector{Float64}; BKT::Vector{Float64}
    DLX::Vector{Float64}; DLXR::Vector{Float64}
    SLOPEC::Vector{Float64}; SINAC::Vector{Float64}

    H::Matrix{Float64}; H1::Matrix{Float64}; H2::Matrix{Float64}
    BH1::Matrix{Float64}; BH2::Matrix{Float64}
    BHR1::Matrix{Float64}; BHR2::Matrix{Float64}
    AVHR::Matrix{Float64}; BHRATIO::Matrix{Float64}

    B::Matrix{Float64}; BI::Matrix{Float64}; BB::Matrix{Float64}
    BH::Matrix{Float64}; BHR::Matrix{Float64}; BR::Matrix{Float64}
    EL::Matrix{Float64}
    AVH1::Matrix{Float64}; AVH2::Matrix{Float64}; BNEW::Matrix{Float64}

    DEPTHB::Matrix{Float64}; DEPTHM::Matrix{Float64}
    FETCHU::Matrix{Float64}; FETCHD::Matrix{Float64}

    Z::Vector{Float64}; ELWS::Vector{Float64}; SELWS::Vector{Float64}
    BCONSTRICTION::Vector{Float64}

    HTMP1::Float64; HTMP2::Float64   # scratch temporaries, kept for parity

    CONSTRICTION::Matrix{Bool}

    # --- grid/waterbody definition block (w2_con.csv) ---
    # NOTE: LAT/LONGIT are declared in MODULE SURFHE in the Fortran source, not
    # GEOMC -- grouped here anyway because they're read in the same w2_con.csv
    # block as the rest of grid definition and consumed immediately by geometry
    # init (solar angle calcs need them, but so does just plotting the grid).
    LAT::Vector{Float64}; LONGIT::Vector{Float64}     # per waterbody (NWB)
    ELBOT::Vector{Float64}                            # bottom elevation, per waterbody (NWB)
    NL::Vector{Int}                                   # number of layers, per branch (NBR)
end

function W2Geometry()
    W2Geometry(
        Int[], Int[], Int[], Int[],
        Float64[], Float64[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[], Float64[],
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0),
        Float64[], Float64[], Float64[],
        Float64[],
        0.0, 0.0,
        zeros(Bool, 0, 0),
        Float64[], Float64[],
        Float64[],
        Int[],
    )
end

# ------------------------------------------------------------------------------
# MODULE NAMESC
# ------------------------------------------------------------------------------
"""
    W2Names

Mirrors `MODULE NAMESC` -- constituent/output naming and formatting tables.
These are populated directly from the w2_con.csv "CST" (constituent state
variable) block once the CSV reader exists.
"""
mutable struct W2Names
    LNAME::Vector{Int}
    CUNIT::Vector{String}; CUNIT2::Vector{String}; CUNIT3::Vector{String}
    CNAME2::Vector{String}; CDNAME2::Vector{String}
    FMTH::Vector{String}; FMTC::Vector{String}; FMTCD::Vector{String}
    CNAME1::Vector{String}; CDNAME1::Vector{String}
    CNAME::Vector{String}; CNAME3::Vector{String}
    CDNAME::Vector{String}; CDNAME3::Vector{String}; HNAME::Vector{String}
    TITLE::Vector{String}
    CONV::Matrix{String}
end

function W2Names()
    W2Names(
        Int[],
        String[], String[], String[],
        String[], String[],
        String[], String[], String[],
        String[], String[],
        String[], String[],
        String[], String[], String[],
        String[],
        fill("", 0, 0),
    )
end

# ------------------------------------------------------------------------------
# Time / stability control (w2_con.csv "Time control cards" + adjacent block)
# ------------------------------------------------------------------------------
"""
    W2TimeControl

Read directly after the dimension block in w2_con.csv, before grid definition.
Fortran source for these fields is split across GLOBAL (DLT family) and RSTART
-- grouped here because they're read together as one control-file section and
nothing downstream needs them split by original module.
"""
mutable struct W2TimeControl
    TMSTRT::Float64; TMEND::Float64; YEAR::Int

    NDLT::Int; DLTMIN::Float64; DLTINTER::Bool   # ON/OFF -> Bool
    DLTD::Vector{Float64}    # length NDLT -- julian day breakpoints for timestep changes
    DLTMAX::Vector{Float64}  # length NDLT -- max timestep at each breakpoint
    DLTF::Vector{Float64}    # length NDLT -- fraction-of-max-stable-step at each breakpoint

    VISC::Vector{Bool}; CELC::Vector{Bool}; DLTADD::Vector{Bool}   # per waterbody (NWB)
end

function W2TimeControl()
    W2TimeControl(0.0, 0.0, 0,
        0, 0.0, false,
        Float64[], Float64[], Float64[],
        Bool[], Bool[], Bool[],
    )
end

end # module W2Core