# ==============================================================================
# Plotting/LongitudinalProfile.jl
#
# A standalone sanity-check plot: the side-view (longitudinal profile) of the
# waterbody -- water surface and channel bottom vs. distance along each
# branch. This is purely a debugging/validation tool for the IO module, not
# part of the simulation core -- kept in its own file/module on purpose so it
# can depend on a plotting package (Plots.jl) without that dependency
# touching Core, IO, or the solvers.
#
# SCOPE NOTE: this is a deliberately rough first cut, not a faithful
# longitudinal-profile renderer. A geometrically exact version needs KB(I)
# (bottom-active-layer index per segment) and EL(K,I) (layer elevations),
# both of which come from init-geom.F90's logic -- not yet ported. Until
# then, this approximates each segment's bottom elevation as:
#     ELBOT(waterbody) + (number of inactive bottom layers, i.e. zero-width
#     rows in B(:,I), counted from the bottom) -- a stand-in, not the real
#     stepped channel bottom W2 actually computes.
# Good enough to catch gross errors (wrong segment count, wrong branch
# boundaries, a bathymetry file read offset by one row) -- not good enough
# to be the figure in a report.
#
# LAYOUT: single row, one axes for the whole grid -- x = segment index
# (1:IMX), not per-branch cumulative distance. Branches are plotted as
# separate line segments (with a NaN gap between them so disjoint branches
# don't get connected by a straight line across the boundary/inactive
# segments between them), sharing one legend for the whole figure rather
# than one legend box per branch.
# ==============================================================================

module LongitudinalProfile

using Plots

export plot_longitudinal_profile

"""
    plot_longitudinal_profile(g, geom, bathy_by_wb)

`g`, `geom` are the structs from InputReader.read_control_file.
`bathy_by_wb` is a Dict{Int, NamedTuple} mapping waterbody index -> the
named tuple returned by BathymetryReader.read_bathymetry! for that waterbody.

Single-row plot: x = segment index (1:IMX), y = elevation (m). Water surface
as a line, bottom (approx.) as a line with a filled area below it. All
branches share one x-axis and one legend box (each label appears once).
"""
function plot_longitudinal_profile(g, geom, bathy_by_wb::Dict)
    p = plot(;
        xlabel="segment", ylabel="elevation (m)",
        title="Longitudinal Profile (all branches)",
        legend=:topright, size=(1400, 500),
        left_margin=12Plots.mm, right_margin=12Plots.mm,
        top_margin=6Plots.mm, bottom_margin=10Plots.mm,
        xlims=(0, g.IMX + 1),
    )

    labeled = false
    for jw in 1:g.NWB
        haskey(bathy_by_wb, jw) || continue
        bathy = bathy_by_wb[jw]

        for jb in g.BS[jw]:g.BE[jw]
            i_first, i_last = g.US[jb], g.DS[jb]
            seg = i_first:i_last
            x = collect(seg)

            top = geom.ELWS[seg]

            # Rough bottom estimate: count active (nonzero-width) layers per
            # segment from B(:,I), subtract that depth from the water surface.
            # See module docstring -- this is the part to replace once
            # KB(I)/EL(K,I) exist.
            active_layers = [count(>(0), geom.B[:, i]) for i in seg]
            avg_layer_height = bathy.H[1:g.KMX] |> h -> sum(h) / max(count(>(0), h), 1)
            bottom = top .- active_layers .* avg_layer_height

            plot!(p, x, top;
                label=(labeled ? "" : "water surface"), linewidth=2, color=:steelblue)
            plot!(p, x, bottom;
                label=(labeled ? "" : "bottom (approx.)"), linewidth=1.5, color=:saddlebrown)
            plot!(p, x, bottom; fillrange=top, fillalpha=0.15, color=:steelblue, label="")
            labeled = true
        end
    end

    labeled || error("plot_longitudinal_profile: no branches plotted -- check bathy_by_wb has an entry for at least one waterbody in 1:g.NWB")

    return p
end

end # module LongitudinalProfile
