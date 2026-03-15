#!/usr/bin/env bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "=== Running benchmark with type-instability fixes (new) ==="
julia --project example/benchmark_ngb.jl
cp example/benchmark_results.csv /tmp/benchmark_new.csv

echo ""
echo "=== Stashing src changes (reverting to baseline) ==="
git stash push -- src/balance.jl src/gravity_treewalk.jl src/ngb.jl src/particle_data.jl src/tree.jl

echo ""
echo "=== Running benchmark without fixes (old/baseline) ==="
julia --project example/benchmark_ngb.jl
cp example/benchmark_results.csv /tmp/benchmark_old.csv

echo ""
echo "=== Restoring type-instability fixes ==="
git stash pop

echo ""
echo "=== Comparison: old vs new ==="
julia --project -e "
using Printf, DelimitedFiles

old = readdlm(\"/tmp/benchmark_old.csv\", ',', skipstart=1)
new = readdlm(\"/tmp/benchmark_new.csv\", ',', skipstart=1)

@printf(\"%-10s | %10s %10s %8s | %10s %10s %8s | %10s %10s %8s\n\",
    \"Npart\", \"build_old\", \"build_new\", \"speedup\",
    \"gather_old\", \"gather_new\", \"speedup\",
    \"scat_old\", \"scat_new\", \"speedup\")
println(repeat(\"-\", 105))
for i in 1:size(old,1)
    @printf(\"%-10d | %8.2f ms %8.2f ms %7.2fx | %8.2f μs %8.2f μs %7.2fx | %8.2f μs %8.2f μs %7.2fx\n\",
        Int(old[i,1]),
        old[i,2], new[i,2], old[i,2]/new[i,2],
        old[i,3], new[i,3], old[i,3]/new[i,3],
        old[i,4], new[i,4], old[i,4]/new[i,4])
end
"
