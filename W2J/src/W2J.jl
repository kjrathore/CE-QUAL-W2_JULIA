"""
    W2J

CE-QUAL-W2 in Julia. A from-scratch reimplementation of CE-QUAL-W2's 2D
laterally-averaged hydrodynamic and water quality model, targeting:

  1. Parallel processing  — multi-core speedup via per-cell kinetics and
     per-segment tridiagonal solves (confirmed embarrassingly parallel).
  2. Differentiable programming — gradient-based calibration via Enzyme.jl.
  3. Architectural flexibility — water quality kinetics decoupled from the
     hydrodynamics/grid layer, so it does not assume any particular spatial
     layout and could plug into a different transport engine later.

See README.md for architecture overview and the project research notes for
the full Fortran call-graph audit and ENTRY-decomposition findings that this
design is based on.

Status: IO/geometry foundation implemented and validated against real Detroit
Reservoir data (see W2J_README.md Progress Tracker). Hydrodynamic core (free-
surface/momentum solve, transport, turbulence closure) not yet implemented;
`density()` (Hydrodynamics/Density.jl) is the first piece of it.
"""
module W2J

include("Core/Architecture.jl")
include("Core/State.jl")
include("Core/Parallel.jl")
include("Core/Grid.jl")
include("Core/InitGeometry.jl")

include("Solvers/Tridiagonal.jl")

include("Hydrodynamics/Density.jl")
include("Hydrodynamics/FreeSurface.jl")
include("Hydrodynamics/Waterbody.jl")
include("Hydrodynamics/Transport.jl")
include("Hydrodynamics/Turbulence.jl")
include("Hydrodynamics/Structures.jl")

include("WaterQuality/RateMultipliers.jl")
include("WaterQuality/Kinetics.jl")

include("IO/InputReader.jl")
include("IO/BathymetryReader.jl")
include("IO/OutputWriter.jl")

include("Plotting/LongitudinalProfile.jl")

include("Simulation.jl")

end # module W2J
