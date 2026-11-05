#!/bin/bash
# Copy the summary tables from ANAEROBIC_RESULTS (the working folder, with the raw
# samples) into results/ of each comment folder. Raw samples and per-model
# intermediate files stay in the working folder.
# Usage, from analysis/anaerobic: ANAEROBIC_RESULTS=<folder> bash paper_figures/collect_results.sh
set -euo pipefail
SRC=${ANAEROBIC_RESULTS:?set ANAEROBIC_RESULTS}
put() {   # put <folder> <path in SRC>...
    local dir=$1; shift
    mkdir -p "$dir/results"
    for p in "$@"; do rm -rf "$dir/results/$p"; mkdir -p "$dir/results/$(dirname "$p")"; cp -r "$SRC/$p" "$dir/results/$p"; done
}
put r1_major3_attribution repro function_ablation function_pairs sectionB curation_ablation.tsv curation_pairs.tsv biomassDR.txt biomassTermVariants.txt
put r1_major5_degeneracy fva analysis_opt0.99999 analysis_sensitivity chains
put r1_major6_gcy1 gcy1_thermodynamics.tsv
put r2_major1_deletions deletion_comparison.tsv ko_panel.tsv adh3_gpd2.tsv engineered_panel.tsv gpd_viability.tsv ko_screen frdKnockouts.txt gpdSinkChain.txt
put r2_major2_fig3 flux_comparison.tsv
put r2_major3_protons protonBudget.txt
put paper_figures supplementary figures
rm -f r1_major3_attribution/results/curation_ablation_*.tsv
