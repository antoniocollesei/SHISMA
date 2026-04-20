#!/bin/bash
# Usage: ./wrapper.sh --data <time_series_file> --ppi <ppi_file> --ct <cell_type> --out <output_name> --nperm <num_permutations> --cores <num_cores>

# Set default values
default_num_permutations=100
default_num_cores=1
default_penalization_method="fdr"
default_penalization_threshold=0.05

run() {
    "$@" > /dev/null 2>&1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --data) time_series_file="$2"; shift 2 ;;
        --ppi) ppi_file="$2"; shift 2 ;;
        --ct) cell_type="$2"; shift 2 ;;
        --out) output_name="$2"; shift 2 ;;
        --mht) penalization_method="$2"; shift 2 ;;
        --alpha) penalization_threshold="$2"; shift 2 ;;
        --nperm) num_permutations="$2"; shift 2 ;;
        --cores) num_cores="$2"; shift 2 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
done

num_permutations="${num_permutations:-$default_num_permutations}"
num_cores="${num_cores:-$default_num_cores}"
penalization_method="${penalization_method:-$default_penalization_method}"
penalization_threshold="${penalization_threshold:-$default_penalization_threshold}"

if [[ -z "$time_series_file" || -z "$ppi_file" || -z "$cell_type" || -z "$output_name" ]]; then
    echo "Usage: $0 --data <time_series_file> --ppi <ppi_file> --ct <cell_type> --out <output_name> [--mht <penalization_method>] [--alpha <penalization_threshold>] [--nperm <num_permutations>] [--cores <num_cores>]"
    exit 1
fi

echo "Analysis of $time_series_file started. PPI being used is $ppi_file"
start_time=$(date +%s)

output_folder=analysis_"$output_name"
data=$PWD/"$output_folder"/shapelet_data
intermediate=$PWD/"$output_folder"/shapelet_intermediate
results=$PWD/"$output_folder"/shapelet_results

run mkdir -p "$output_folder"
run mkdir -p "$data"
run mkdir -p "$intermediate"
run mkdir -p "$results"

echo "Generating gene-shape scores..."

run python py/make_penalized_shap.py \
        --time_series "$time_series_file" \
        --ppi "$ppi_file" \
        --cell_type "$cell_type" \
        --output_folder "$output_folder" \
        --penalization_method "$penalization_method" \
        --penalization_threshold "$penalization_threshold" \
        --num_permutations "$num_permutations" \
        --num_cores "$num_cores"

echo "Completed."

network=$(basename "$ppi_file" | sed 's/\.[^.]*$//')
run mkdir -p "$intermediate/$network"

score_files=($(ls "$data"/scores_*.tsv 2>/dev/null | xargs -n 1 basename | sed 's/.tsv//'))

for score_file in "${score_files[@]}"
do
        run mkdir -p "$intermediate/${network}_${score_file}"
done

echo "Constructing similarity matrices for network diffusion..."

run python py/src-hierarchical-hotnet/construct_similarity_matrix.py \
        -i   "$PWD/$output_folder/${network}_edge_list.tsv" \
        -o   "$intermediate/$network/similarity_matrix.h5" \
        -bof "$intermediate/$network/beta.txt"

echo "Completed."

echo "Permuting networks..."

run cp "$PWD/$output_folder/${network}_index_gene.tsv" "$intermediate/$network/index_gene_0.tsv"
run cp "$PWD/$output_folder/${network}_edge_list.tsv" "$intermediate/$network/edge_list_0.tsv"

run parallel -u -j "$num_cores" \
        python py/src-hierarchical-hotnet/permute_network.py \
                -i "$intermediate/$network/edge_list_0.tsv" \
                -s {} \
                -c \
                -o "$intermediate/$network/edge_list_{}.tsv" \
        ::: $(seq 1 4)

run parallel -u -j "$num_cores" \
        python py/src-hierarchical-hotnet/permute_network.py \
                -i "$intermediate/$network/edge_list_0.tsv" \
                -s {} \
                -o "$intermediate/$network/edge_list_{}.tsv" \
        ::: $(seq 5 8)

echo "Completed."

echo "Permuting scores..."

for score_file in "${score_files[@]}"
do
        run cp "$data/$score_file.tsv" "$intermediate/${network}_${score_file}/scores_0.tsv"

        run python py/src-hierarchical-hotnet/find_permutation_bins.py \
                -gsf "$intermediate/${network}_${score_file}/scores_0.tsv" \
                -igf "$PWD/$output_folder/${network}_index_gene.tsv" \
                -elf "$PWD/$output_folder/${network}_edge_list.tsv" \
                -ms  1000 \
                -o   "$intermediate/${network}_${score_file}/score_bins.tsv"

        run parallel -u -j "$num_cores" \
                python py/src-hierarchical-hotnet/permute_scores.py \
                        -i  "$intermediate/${network}_${score_file}/scores_0.tsv" \
                        -bf "$intermediate/${network}_${score_file}/score_bins.tsv" \
                        -s  {} \
                        -o  "$intermediate/${network}_${score_file}/scores_{}.tsv" \
                ::: $(seq "$num_permutations")
done

echo "Completed."

echo "Constructing hierarchies..."

for score_file in "${score_files[@]}"
do
        run parallel -u -j "$num_cores" \
                python py/src-hierarchical-hotnet/construct_hierarchy.py \
                        -smf  "$intermediate/$network/similarity_matrix.h5" \
                        -igf  "$PWD/$output_folder/${network}_index_gene.tsv" \
                        -gsf  "$intermediate/${network}_${score_file}/scores_{}.tsv" \
                        -helf "$intermediate/${network}_${score_file}/hierarchy_edge_list_{}.tsv" \
                        -higf "$intermediate/${network}_${score_file}/hierarchy_index_gene_{}.tsv" \
                ::: $(seq 0 "$num_permutations")
done

echo "Processing hierarchies..."

for score_file in "${score_files[@]}"
do
        target_file="$intermediate/${network}_${score_file}/hierarchy_edge_list_0.tsv"

        if [[ -f "$target_file" ]]; then
                run python py/src-hierarchical-hotnet/process_hierarchies.py \
                        -oelf "$intermediate/${network}_${score_file}/hierarchy_edge_list_0.tsv" \
                        -oigf "$intermediate/${network}_${score_file}/hierarchy_index_gene_0.tsv" \
                        -pelf $(for i in $(seq "$num_permutations"); do echo "$intermediate/${network}_${score_file}/hierarchy_edge_list_${i}.tsv"; done) \
                        -pigf $(for i in $(seq "$num_permutations"); do echo "$intermediate/${network}_${score_file}/hierarchy_index_gene_${i}.tsv"; done) \
                        -lsb  1 \
                        -cf   "$results/clusters_${network}_${score_file}.tsv" \
                        -pl   "$network" "$score_file" \
                        -pf   "$results/sizes_${network}_${score_file}.pdf" \
                        -nc   "$num_cores"
        else
                echo "Skipping $score_file because hierarchy not being evaluated for this shape."
        fi
done

echo "Completed."

echo "Parsing results..."

run python py/parse_hotnet_results.py --time_series "$time_series_file" --cell_type "$cell_type" --results_dir "$results" --output_folder "$output_folder"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "SHISMA's execution completed in $elapsed_time seconds."
