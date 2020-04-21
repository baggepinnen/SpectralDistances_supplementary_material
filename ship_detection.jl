using Plots, LPVSpectral, DSP, SpectralDistances, ControlSystems, TotalLeastSquares, Random, WAV, Optim
using Base.Iterators: product
using IterTools: subsets

Continuous = SpectralDistances.Continuous
default(grid = false)
cd(@__DIR__)
y = [wavread("ships/y$(i).wav")[1][:] for i in 1:3]

fitmethod = SpectralDistances.TLS(na = 4)
models = fitmethod.(y)

## Barycenter ===================================================
options = (
    solver     = sinkhorn_log!,
    tol        = 1e-8,
    iters      = 1_000_000,
    γ          = 0.0,
    uniform    = true,
    inneriters = 500_000,
    innertol   = 1e-6,
)
distance = OptimalTransportRootDistance(
    domain     = Continuous(),
    p          = 2,
    β          = 0.01,
    weight     = unitweight,
)

Random.seed!(1)
bc = barycenter(distance, models; options...)
w  = exp10.(LinRange(-0.5, 0.5, 550))
plot(legend = :bottomleft)
bodeplot!.( models, Ref(w),
    plotphase = false,
    lab       = "Input",
    color     = :blue,
    linestyle = :auto,
)
bodeplot!( bc, w,
    plotphase = false,
    lab       = "Barycenter",
    xscale    = :identity,
    title     = "",
    color     = :orange,
    xticks    = false,
    yticks    = false,
    grid      = false,
)



## Barycentric coordinates ==============================

# Create mix signal
ymix     = sum(y[1:2])
modelmix = TLS(na = 8)(ymix)

fig_bode = bodeplot!( modelmix, w,
    plotphase = false,
    lab       = "Mix signal",
    xscale    = :identity,
    xticks    = false,
    yticks    = false,
    grid      = false,
    color     = :cyan,
    xlabel    = "Frequency",
    title     = "Rational DEMON Spectra",
    ylabel    = "",
)

## Calculate barycentric coordinates
using Optim
optimopts = Optim.Options(
    store_trace       = true,
    show_trace        = true,
    show_every        = 1,
    iterations        = 20,
    allow_f_increases = false,
    time_limit        = 100,
    x_tol             = 1e-5,
    f_tol             = 1e-6,
    g_tol             = 1e-6,
)
options = (
    tol        = 1e-7,
    iters      = 200_000,
    γ          = 1,
    uniform    = true,
    inneriters = 200_000,
    innertol   = 1e-6,
)
meth = ParticleSwarm()
λ = barycentric_coordinates( distance, models, modelmix, meth;
    options = optimopts,
    plot    = nothing,
    robust  = false,
    uniform = false,
    options...,
)
println("λ: ", round.(λ, digits = 5))
fig_bar = bar(λ,
    title  = "Barycentric Coordinates \$\\lambda\$",
    legend = false,
    xticks = 1:3,
    xlabel = "Ship Number",
)


## Consider all poles for transport
# This is the code to reproduce Figure 11

ymix     = sum(y[1:2]) + sin.(2π * 0.25 .* (1:size(y[1], 1))) # Add noise component to mix signal
modelmix = fitmodel(SpectralDistances.TLS(na = 10, λ = 1e-1), ymix)

distance = OptimalTransportRootDistance(
    domain     = Continuous(),
    p          = 2,
    β          = 0.01,
    weight     = e -> ones(size(e)),
    divergence = KL(0.01),
)
allroots = reduce(vcat, [m.pc for m in models])
xlabs    = ["$i:$j" for i = 1:3 for j = 1:4]
xt       = (1:length(xlabs), xlabs)
Γ        = transport_plan(distance, modelmix.pc, allroots)
heatmap(Γ,
    xtick  = xt,
    color  = :viridis,
    xlabel = "Signal number : Pole number",
    ylabel = "Mixsignal pole number",
)

##
plot(legend = :bottomleft)
bodeplot!.( models, Ref(w),
    plotphase = false,
    lab       = "Input",
    color     = :blue,
    linestyle = :auto,
)
bodeplot!(modelmix, w,
    plotphase = false,
    lab       = "Mix signal",
    xscale    = :identity,
    xticks    = false,
    yticks    = false,
    grid      = false,
    color     = :cyan,
    xlabel    = "Frequency",
    title     = "Rational DEMON Spectra",
    ylabel    = "",
)
