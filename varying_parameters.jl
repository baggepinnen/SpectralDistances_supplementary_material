using SpectralDistances, DSP, Statistics, LinearAlgebra, Distances, Random, LaTeXStrings
const SD     = SpectralDistances
rescale(x)   = x ./ mean(maximum(x))
default(grid = false)

function plotall(generate_signal, X, fvec; order = 4, p = 2, kwargs...)

    loss = ModelDistance(
        LS(na = order),
        RationalOptimalTransportDistance(domain = Continuous(), p = p),
    )
    losses_closedform = map(fvec) do f
        SD.evaluate(loss, X, generate_signal(f))
    end
    plot( fvec, rescale(losses_closedform);
        subplot    = 1,
        lab        = L"W_2",
        legendfont = font(5),
        layout     = 2,
        yscale     = :log10,
        kwargs...,
    )
    plot!( fvec, rescale(losses_closedform);
        legendfont = font(5),
        xscale     = get(kwargs, :xscale, :identity),
        xlabel     = get(kwargs, :xlabel, ""),
        subplot    = 2,
        lab        = "",
        legend     = false,
        kwargs...,
    )

    loss = WelchOptimalTransportDistance(distmat = nothing, args = (128,), p = p)
    losses_pwelch = map(fvec) do f # segfaults with threads
        SD.evaluate(loss, X, generate_signal(f))
    end
    plot!(fvec, clamp.(rescale(losses_pwelch), 1e-5, 1); lab = L"W_2 Welch", l = :dash)
    plot!(fvec, rescale(losses_pwelch); subplot = 2, lab = "", l = :dash)

    loss = ModelDistance(
        LS(na = order),
        EuclideanRootDistance(domain = Continuous(), weight = s1 ∘ residueweight, p = p),
    )
    losses_roots_cont = map(fvec) do f
        SD.evaluate(loss, X, generate_signal(f))
    end
    plot!(fvec, rescale(losses_roots_cont), subplot = 1, lab = "WRD", l = :auto)
    plot!(fvec, rescale(losses_roots_cont), subplot = 2, lab = "", l = :auto)


    loss = ModelDistance(
        LS(na = order),
        SinkhornRootDistance(domain = Continuous(), β = 0.0001, p = p),
    )
    losses_sinkhorn = map(fvec) do f
        SD.evaluate( loss, X,
            generate_signal(f),
            solver = sinkhorn_log!,
            tol    = 1e-9,
            iters  = 500_000,
        )
    end
    plot!(fvec, rescale(losses_sinkhorn), subplot = 1, lab = "OTRD", l = :dash)
    plot!(fvec, rescale(losses_sinkhorn), subplot = 2, lab = "", l = :dash)

    display(current())
end


# # Interpolate frequencies
t               = 0:0.05:1000
fvec            = exp10.(LinRange(-2, -0.05, 100))
generate_signal = f -> sin.(2pi * f .* t)
X               = sin.(2pi * 0.1 .* t)
plotall( generate_signal, X, fvec;
    order  = 2,
    p      = 2,
    xscale = :log10,
    xlabel = "Frequency",
    legend = true,
)


# # Different cutoff lowpass filters
function bp_filter(x, passband)
    responsetype = Bandpass(passband..., fs = 1)
    designmethod = Butterworth(2)
    filt(digitalfilter(responsetype, designmethod), x)
end

Random.seed!(123)
let fvec = exp10.(LinRange(-1.5, log10(0.45), 100)),
    X               = bp_filter(randn(10000), (1e-3, 0.1)),
    x0              = randn(10000),
    generate_signal = f -> bp_filter(x0, (1e-3, f))

    plotall( generate_signal, X, fvec,
        order  = 6,
        p      = 2,
        xlabel = "Cutoff Frequency",
        xscale = :log10,
    )
end


# # Different cutoff highpass filters
Random.seed!(123)
let fvec = exp10.(LinRange(-2.5, log10(0.4), 100)),
    X               = bp_filter(randn(10000), (0.1, 0.45)),
    x0              = randn(10000),
    generate_signal = f -> bp_filter(x0, (f, 0.45))

    @time plotall( generate_signal, X, fvec,
        order     = 6,
        p         = 2,
        xlabel    = "Cutoff Frequency",
        xscale    = :log10,
        fillalpha = 0.1,
    )
end
