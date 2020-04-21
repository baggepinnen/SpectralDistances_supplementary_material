# # Illustration of cumulative function
using ControlSystems, SpectralDistances, LaTeXStrings
default(grid = false)

color_palette = [:blue, :orange, :magenta, :cyan]

G1   = tf(1, [1, 0.12, 1]) * tf(1, [1, 0.1, 0.1])
G2   = tf(1, [1, 0.12, 2]) * tf(1, [1, 0.1, 0.4])
a1   = denvec(G1)[]
a2   = denvec(G2)[]
n    = length(a1)

f1c  = w -> abs2(1 / sum(j -> a1[j] * (im * w)^(n - j), 1:n))
f2c  = w -> abs2(1 / sum(j -> a2[j] * (im * w)^(n - j), 1:n))
sol1 = SpectralDistances.c∫(f1c, 0, 2π)
sol2 = SpectralDistances.c∫(f2c, 0, 2π)
plot(
    (sol1.t .+ sol1.t[2]) .* 2π,
    sqrt.(sol1 ./ sol1[end]),
    label     = L"F_1(w)",
    fillrange = sqrt.(sol2(sol1.t) ./ sol2[end]),
    fill      = (0.6, color_palette[3]),
    l         = (2, color_palette[1]),
)
plot!(
    (sol2.t .+ sol2.t[2]) .* 2π,
    sqrt.(sol2(sol2.t) ./ sol2[end]),
    label  = L"F_2(w)",
    l      = (2, color_palette[2]),
    xscale = :log10,
    legend = false,
    grid   = false,
    xlabel = L"\omega",
    xlims  = (1e-2, 2pi),
)

##

bodeplot(
    [G1, G2],
    exp10.(LinRange(-1.5, 1, 200)),
    legend    = false,
    grid      = false,
    title     = "",
    linecolor = [color_palette[1] color_palette[2]],
    l         = (2,),
    plotphase = false,
)

##

pzmap( [G1, G2],
    legend      = false,
    grid        = false,
    title       = "",
    markercolor = [color_palette[1] color_palette[2]],
    color       = [color_palette[1] color_palette[2]],
    m           = (2, :c),
    xlims       = (-0.5, 0.5),
    xlabel      = "Re",
    ylabel      = "Im",
)
vline!([0], l = (:black, :dash))
hline!([0], l = (:black, :dash))
