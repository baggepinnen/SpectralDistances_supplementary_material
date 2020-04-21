using SpectralDistances, ControlSystems, Distances, Plots, LaTeXStrings, Random
using SpectralDistances: domain, interpolator
Random.seed!(2) # Choose a seed that produces visually pleasing spectra
n      = 4
r1     = complex.(-0.01 .+ 0.001 * randn(3), 2 * randn(3))
r1     = ContinuousRoots([r1; conj.(r1)])

r2     = complex.(-0.01 .+ 0.001 * randn(3), 2 * randn(3))
r2     = ContinuousRoots([r2; conj.(r2)])
r1, r2 = normalize_energy.((r1, r2))
A1     = AR(r1)
A2     = AR(r2)

## Interpolation of spectra
fig1   = plot()
t      = 0.1
dist   = ClosedFormSpectralDistance(domain = Continuous(), p = 2, interval = (0.0, exp10(1.6)))
interp = SpectralDistances.interpolator(dist, A1, A2)
w      = exp10.(LinRange(-1.5, 1, 300))
for t in LinRange(0, 1, 7)
    Φ = clamp.(interp(w, t), 1e-10, 100)
    plot!( w, sqrt.(Φ),
        xscale   = :log10,
        yscale   = :log10,
        line_z   = t,
        lab      = "",
        xlabel   = "",
        title    = L"W_2",
        ylims    = (1e-3, 1e1),
        colorbar = false,
        l        = (1,),
        c        = :viridis,
    )
end

fig2 = plot()
w    = exp10.(LinRange(-1.5, 1, 300))
Φ1   = bode(tf(A1), w)[1][:]
Φ2   = bode(tf(A2), w)[1][:]
for t in LinRange(0, 1, 7)
    plot!(
        w,
        (1 - t) .* Φ1 .+ t .* Φ2,
        xscale   = :log10,
        yscale   = :log10,
        line_z   = t,
        lab      = "",
        xlabel   = "",
        title    = L"L_2",
        ylims    = (1e-3, 1e1),
        colorbar = false,
        l        = (1,),
        c        = :viridis,
    )
end

##
rdist  = EuclideanRootDistance(domain = Continuous(), p = 2, weight = e -> ones(length(e)))
interp = SpectralDistances.interpolator(rdist, A1, A2, normalize = false)
fig3   = plot()
w      = exp10.(LinRange(-1.5, 1, 300))
for t in LinRange(0, 1, 7)
    Φ = interp(w, t)
    plot!( fig3, w, sqrt.(Φ),
        xscale   = :log10,
        yscale   = :log10,
        line_z   = t,
        lab      = "",
        xlabel   = "",
        title    = L"RD",
        ylims    = (1e-3, 1e1),
        colorbar = false,
        l        = (1,),
        c        = :viridis,
    )
end

##
rdist  = EuclideanRootDistance(domain = Continuous(), p = 2, weight = residueweight)
interp = SpectralDistances.interpolator(rdist, A1, A2, normalize = false)
fig4   = plot()
w      = exp10.(LinRange(-1.5, 1, 300))
for t in LinRange(0, 1, 7)
    Φ = interp(w, t)
    plot!( w, sqrt.(Φ),
        xscale   = :log10,
        yscale   = :log10,
        line_z   = t,
        lab      = "",
        xlabel   = "",
        title    = L"WRD",
        ylims    = (1e-3, 1e1),
        colorbar = false,
        l        = (1,),
        c        = :viridis,
    )
end


fig = plot(fig1, fig2, fig3, fig4, layout = (2, 2))
display(fig)


## Pole interpolation
function pole_interpolator(d::EuclideanRootDistance, A1, A2)
    p = d.p
    @assert p == 2 "Interpolation only supported for p=2, you have p=$p"
    RT = domain(d) isa Continuous ? ContinuousRoots : DiscreteRoots
    e1, e2 = roots(domain(d), A1), roots(domain(d), A2)
    I1, I2 = d.assignment(e1, e2)
    e1, e2 = RT(e1[I1]), RT(e2[I2])
    w1, w2 = d.weight(e1), d.weight(e2)
    function (t)
        RT(((1 - t) * w1 .* e1 + t * w2 .* e2) ./ ((1 - t) .* w1 .+ t .* w2))
    end
end

pdist      = EuclideanRootDistance(domain = Continuous(), p = 2, weight = unitweight)
poleinterp = pole_interpolator(pdist, A1, A2)
figpoles   = plot()
for t in LinRange(0, 1, 30)
    plot!(
        figpoles,
        poleinterp(t),
        marker_z = t,
        lab      = "",
        xlabel   = "",
        title    = "",
        m        = (5,),
        circle   = t == 0,
        c        = :viridis,
        xlims    = (-0.05, 0),
        colorbar = false,
        xlabel   = "Re",
        ylabel   = "Im",
        markerstrokealpha = 0,
    )
end
hline!([0], l = (:black, :dash), primary = false)
vline!([0], l = (:black, :dash), primary = false)
