## Generate random clusters
using SpectralDistances, Random, DSP
Random.seed!(1)
Discrete = SpectralDistances.Discrete
# Generate signals
N_clusters = 50
N_signals  = 50
signals = map(1:N_clusters) do i
    passband = sort(0.5 * rand(2))
    signals  = map(1:N_signals) do s
        SpectralDistances.bp_filter(randn(500), passband)
    end

end
assi      = repeat((1:N_clusters)', N_signals)[:]
allsounds = reduce(vcat, signals)

# Fit models
na        = 4
fitmethod = TLS(na = na)
models    = map(allsounds) do sound
    fitmodel(fitmethod, sound)
end


embeddings = ContinuousRoots.(move_real_poles.(roots.(Ref(Discrete()), models), 1e-2))
X0         = reduce(hcat, embeddings)
X          = Float64.([real(X0[1:end÷2, :]); imag(X0[1:end÷2, :])])
X          = SpectralDistances.m1(X, 1)

# Fit welch periodograms
periodograms = welch_pgram.(allsounds, 128)
pembeddings  = power.(periodograms)
pX           = reduce(hcat, pembeddings)
pX           = SpectralDistances.m1(pX, 1)

## Plot
scatter(twoD(X'),
    marker_z          = assi,
    m                 = (5,),
    colorbar          = false,
    markerstrokealpha = 0,
    c                 = :viridis,
    layout            = 2,
    title             = "Root embedding",
    legend            = false,
    axis              = false,
    titlefont         = font("times", 25),
)
scatter!(twoD(pX'),
    marker_z          = assi,
    m                 = (5,),
    colorbar          = false,
    markerstrokealpha = 0,
    c                 = :viridis,
    subplot           = 2,
    title             = "Welch embedding",
    legend            = false,
    axis              = false,
    titlefont         = font("times", 25),
)
