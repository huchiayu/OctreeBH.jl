using OctreeBH
using StaticArrays
using BenchmarkTools
using Printf
using Statistics
using Random
using UnicodePlots
using PyPlot

const N = 3
const T = Float64
const BOXSIZE = 1.0
const Nngb_target = 32

function format_npart(n)
    s = string(n)
    parts = String[]
    while length(s) > 3
        pushfirst!(parts, s[end-2:end])
        s = s[1:end-3]
    end
    pushfirst!(parts, s)
    return join(parts, "_")
end

Random.seed!(42)

Npart_list = [1_000, 3_000, 10_000, 30_000, 100_000]
Nqueries = 100

# Storage for scaling plot
build_times_ms   = Float64[]
gather_times_us  = Float64[]
scatter_times_us = Float64[]
mean_nngbs       = Float64[]

println("Npart      | build (ms) | gather (μs/query) | scatter (μs/query) | mean Nngb")
println("           |            |     [min time]    |     [min time]     |          ")
println("-----------|------------|-------------------|--------------------|----------")

for Npart in Npart_list
    hsml0 = (Nngb_target * 3 / (4π * Npart))^(1/3) * BOXSIZE

    positions = [SVector{N,T}(rand(), rand(), rand()) for _ in 1:Npart]
    part = [Data{N,T}(positions[i], i, hsml0, 1.0) for i in 1:Npart]

    center         = SVector{N,T}(fill(BOXSIZE / 2, N))
    topnode_length = SVector{N,T}(fill(BOXSIZE, N))
    boxsizes       = SVector{N,T}(fill(BOXSIZE, N))
    query_points   = [SVector{N,T}(rand(), rand(), rand()) for _ in 1:Nqueries]

    print("  Npart=$(format_npart(Npart)): warmup..."); flush(stdout)
    tree = buildtree(part, center, topnode_length)
    for q in query_points
        get_gather_ngb_tree(q, hsml0, tree, boxsizes)
        get_scatter_ngb_tree(q, tree, boxsizes)
    end

    print(" build..."); flush(stdout)
    # Benchmark tree build (single @belapsed — only times one thing)
    build_time = @belapsed buildtree($part, $center, $topnode_length)

    # Rebuild tree for searches
    tree = buildtree(part, center, topnode_length)

    print(" gather..."); flush(stdout)
    gather_min = @belapsed for q in $query_points
        get_gather_ngb_tree(q, $hsml0, $tree, $boxsizes)
    end

    print(" scatter..."); flush(stdout)
    scatter_min = @belapsed for q in $query_points
        get_scatter_ngb_tree(q, $tree, $boxsizes)
    end
    println(" done.")

    gather_median  = gather_min / Nqueries
    scatter_median = scatter_min / Nqueries

    mean_nngb = mean(length(get_gather_ngb_tree(q, hsml0, tree, boxsizes)) for q in query_points)

    push!(build_times_ms,   build_time * 1e3)
    push!(gather_times_us,  gather_median * 1e6)
    push!(scatter_times_us, scatter_median * 1e6)
    push!(mean_nngbs,       mean_nngb)

    @printf("%9s  |    %7.1f  |            %6.1f  |              %6.1f  |    %5.1f\n",
        format_npart(Npart),
        build_time * 1e3,
        gather_median * 1e6,
        scatter_median * 1e6,
        mean_nngb)
end

# Scaling plots
println()
p1 = lineplot(Npart_list, build_times_ms;
    xlabel="Npart", ylabel="build time (ms)",
    title="Tree build time scaling",
    xscale=:log10, yscale=:log10,
    name="build")
println(p1)

p2 = lineplot(Npart_list, gather_times_us;
    xlabel="Npart", ylabel="time per query (μs)",
    title="Neighbor search time scaling",
    xscale=:log10, yscale=:log10,
    name="gather")
lineplot!(p2, Npart_list, scatter_times_us; name="scatter")
println(p2)

# Save PNG with PyPlot
fig, (ax1, ax2) = subplots(1, 2, figsize=(10, 4))

Ns = Float64.(Npart_list)

# O(N log N) reference, normalized to first build point
nlogn_ref = @. Ns * log(Ns)
nlogn_ref .*= build_times_ms[1] / nlogn_ref[1]

# O(log N) reference, normalized to first gather point
logn_ref = log.(Ns)
logn_ref .*= gather_times_us[1] / logn_ref[1]

ax1.loglog(Npart_list, build_times_ms, "o-", label="measured")
ax1.loglog(Npart_list, nlogn_ref, "k--", alpha=0.5, label="O(N log N)")
ax1.set_xlabel("Npart")
ax1.set_ylabel("build time (ms)")
ax1.set_title("Tree build time scaling")
ax1.legend()
ax1.grid(true, which="both", ls="--", alpha=0.4)

ax2.loglog(Npart_list, gather_times_us,  "o-", label="gather")
ax2.loglog(Npart_list, scatter_times_us, "s-", label="scatter")
ax2.loglog(Npart_list, logn_ref, "k--", alpha=0.5, label="O(log N)")
ax2.set_xlabel("Npart")
ax2.set_ylabel("time per query (μs)")
ax2.set_title("Neighbor search time scaling")
ax2.legend()
ax2.grid(true, which="both", ls="--", alpha=0.4)

tight_layout()
PyPlot.savefig("example/benchmark_ngb_scaling.png", dpi=150)
close(fig)
println("\nScaling plot saved to example/benchmark_ngb_scaling.png")

open("example/benchmark_results.csv", "w") do io
    println(io, "npart,build_ms,gather_us,scatter_us,mean_nngb")
    for i in eachindex(Npart_list)
        @printf(io, "%d,%.4f,%.4f,%.4f,%.2f\n",
            Npart_list[i], build_times_ms[i], gather_times_us[i],
            scatter_times_us[i], mean_nngbs[i])
    end
end
println("Results written to example/benchmark_results.csv")
