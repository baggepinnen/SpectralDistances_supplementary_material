cd(@__DIR__)

using WAV, SpectralDistances, AudioClustering, TotalLeastSquares, ControlSystems

path          = "birds/test_padded_30birds/"
files         = readdir(path)
labels0       = getfield.(match.(r"[a-z_]+", files), :match) .|> String
files         = path .* files

ulabels       = unique(labels0)
n_classes_tot = length(ulabels)
fs            = 44100

n_classes     = 30
n_per_class   = 20
labels = sum((labels0 .== reshape(ulabels, 1, :)) .* (1:n_classes_tot)', dims = 2)[:]

classinds = 1:n_classes
n_classes = length(classinds)

# Subsample the dataset to make it smaller
selected_inds = reduce(vcat, [findall(labels .== l)[1:n_per_class] for l in classinds])
selected_files = files[selected_inds]
selected_labels = labels[selected_inds]
selected_birds = getindex.(match.(r"/(\w+)_\d*.wav", selected_files), 1)
selected_birds = replace.(selected_birds, '_' => ' ')
sounds = [wavread(f)[1][:] for f in selected_files]

## Fit rational models (spectra) to the sound signals
fitmethod = TLS(na = 12)
models = SpectralDistances.fitmodel.(fitmethod, sounds)
G = tf.(models)

## Run the clustering algorithm
using Clustering
dist = OptimalTransportRootDistance(domain = Continuous(), β = 0.015, weight = residueweight)
@time clusterresult = SpectralDistances.kbarycenters(
    dist,
    models,
    n_classes,
    seed    = :rand,
    solver  = SpectralDistances.sinkhorn_log!,
    tol     = 1e-4,
    iters   = 2000,
    verbose = true,
    output  = :best,
    uniform = true,
    kiters  = 10,
)

bc, ass = clusterresult.barycenters, clusterresult.assignments

using MLBase, Plots.PlotMeasures
newass, perm = AudioClustering.associate_clusters(selected_labels, ass)
yt = (classinds, [selected_birds[findfirst(selected_labels .== i)] for i in classinds])
@show mean(selected_labels .== newass)

##
cm = confusmat(n_classes, selected_labels, newass)
heatmap(cm ./ sum(cm, dims = 2),
    xlabel         = "Cluster assignment",
    ylabel         = "Best matching class",
    title          = "Confusion Matrix",
    tickfontfamily = "times",
    tickfontsize   = 12,
    margin         = 5mm,
    size           = (1000, 900),
    color          = :viridis,
)
anns = [
    (reverse(ci.I)..., text(val, 12, :gray))
    for (ci, val) in zip(CartesianIndices(cm)[:], vec(cm))
]
annotate!(anns)
yticks!(yt)
xticks!(yt, xrotation = 45)
current()

## Plot bode diagrams of the barycenters
w = exp10.(LinRange(-2, 0.3, 100))
Gbc = tf.(bc)
plot(legend = false, layout = n_classes)
@time for i = 1:n_classes
    Gi = tf(barycenter(
        dist,
        models[selected_labels.==i],
        solver = SpectralDistances.sinkhorn_log!,
    ))
    bodeplot!.(G[selected_labels.==i], Ref(w), plotphase = false, sp = i)
    bodeplot!(Gi, w, plotphase = false, sp = i, l = (:blue))
end
current()
