#!/bin/bash

# Usage: ./wrapper.sh <time_series_file> <ppi_file> <cell_type> <output_folder>

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <time_series_file> <ppi_file> <cell_type>"
    exit 1
fi

start_time=$(date +%s)

data=$PWD/data/shapelet_data
intermediate=$PWD/data/shapelet_intermediate
results=$PWD/data/shapelet_results

mkdir -p $data
mkdir -p $intermediate
mkdir -p $results

# Set the correct Python path from your virtual environment
#python py/make_penalized_shap.py --time_series "$1" --ppi "$2" --cell_type "$3"

num_permutations=100
num_cores=5

# Compile Fortran module.
cd py/py/src-hierarchical-hotnet
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ../..

################################################################################
#   Prepare data.
################################################################################

network=$(basename "$2" | sed 's/\.[^.]*$//')
mkdir -p $intermediate/"$network"

# Identify available score files
score_files=($(ls $data/scores_*.tsv 2>/dev/null | xargs -n 1 basename | sed 's/.tsv//'))

for score_file in "${score_files[@]}"
do
    mkdir -p $intermediate/"$network"_"$score_file"
done

################################################################################
#   Construct similarity matrices.
################################################################################

echo "Construct similarity matrices..."

python py/src-hierarchical-hotnet/construct_similarity_matrix.py \
    -i   $PWD/data/"$network"_edge_list.tsv \
    -o   $intermediate/"$network"/similarity_matrix.h5 \
    -bof $intermediate/"$network"/beta.txt

################################################################################
#   Permute data.
################################################################################

echo "Permuting networks..."

cp $PWD/data/"$network"_index_gene.tsv $intermediate/"$network"/index_gene_0.tsv
cp $PWD/data/"$network"_edge_list.tsv $intermediate/"$network"/edge_list_0.tsv

parallel -u -j $num_cores --bar \
    python py/src-hierarchical-hotnet/permute_network.py \
        -i $intermediate/"$network"/edge_list_0.tsv \
        -s {} \
        -c \
        -o $intermediate/"$network"/edge_list_{}.tsv \
    ::: $(seq 1 4)

parallel -u -j $num_cores --bar \
    python py/src-hierarchical-hotnet/permute_network.py \
        -i $intermediate/"$network"/edge_list_0.tsv \
        -s {} \
        -o $intermediate/"$network"/edge_list_{}.tsv \
    ::: $(seq 5 8)

echo "Permuting scores..."

for score_file in "${score_files[@]}"
do
    cp $data/"$score_file".tsv $intermediate/"$network"_"$score_file"/scores_0.tsv

    python py/src-hierarchical-hotnet/find_permutation_bins.py \
        -gsf $intermediate/"$network"_"$score_file"/scores_0.tsv \
        -igf $data/"$network"_index_gene.tsv \
        -elf $data/"$network"_edge_list.tsv \
        -ms  1000 \
        -o   $intermediate/"$network"_"$score_file"/score_bins.tsv

    parallel -u -j $num_cores --bar \
        python py/src-hierarchical-hotnet/permute_scores.py \
            -i  $intermediate/"$network"_"$score_file"/scores_0.tsv \
            -bf $intermediate/"$network"_"$score_file"/score_bins.tsv \
            -s  {} \
            -o  $intermediate/"$network"_"$score_file"/scores_{}.tsv \
        ::: $(seq $num_permutations)
done

################################################################################
#   Construct hierarchies.
################################################################################

echo "Constructing hierarchies..."

for score_file in "${score_files[@]}"
do
    parallel -u -j $num_cores --bar \
        python py/src-hierarchical-hotnet/construct_hierarchy.py \
            -smf  $intermediate/"$network"/similarity_matrix.h5 \
            -igf  $data/"$network"_index_gene.tsv \
            -gsf  $intermediate/"$network"_"$score_file"/scores_{}.tsv \
            -helf $intermediate/"$network"_"$score_file"/hierarchy_edge_list_{}.tsv \
            -higf $intermediate/"$network"_"$score_file"/hierarchy_index_gene_{}.tsv \
        ::: $(seq 0 $num_permutations)
done

################################################################################
#   Process hierarchies.
################################################################################

echo "Processing hierarchies..."

for score_file in "${score_files[@]}"
do
    python py/src-hierarchical-hotnet/process_hierarchies.py \
        -oelf $intermediate/"$network"_"$score_file"/hierarchy_edge_list_0.tsv \
        -oigf $intermediate/"$network"_"$score_file"/hierarchy_index_gene_0.tsv \
        -pelf $(for i in $(seq $num_permutations); do echo "$intermediate/"$network"_"$score_file"/hierarchy_edge_list_"$i".tsv"; done) \
        -pigf $(for i in $(seq $num_permutations); do echo "$intermediate/"$network"_"$score_file"/hierarchy_index_gene_"$i".tsv"; done) \
        -lsb  1 \
        -cf   $results/clusters_"$network"_"$score_file".tsv \
        -pl   $network $score_file \
        -pf   $results/sizes_"$network"_"$score_file".pdf \
        -nc   $num_cores
done

################################################################################
#   Perform consensus.
################################################################################

echo "Performing consensus..."

python py/src-hierarchical-hotnet/perform_consensus.py \
    -cf  $(for score_file in "${score_files[@]}"; do echo "$results/clusters_"$network"_"$score_file".tsv"; done) \
    -igf $(for score_file in "${score_files[@]}"; do echo "$data/"$network"_index_gene.tsv"; done) \
    -elf $(for score_file in "${score_files[@]}"; do echo "$data/"$network"_edge_list.tsv"; done) \
    -n   $(for score_file in "${score_files[@]}"; do echo "$network"; done) \
    -s   $(for score_file in "${score_files[@]}"; do echo "$score_file"; done) \
    -t   2 \
    -cnf $results/consensus_nodes.tsv \
    -cef $results/consensus_edges.tsv

### Parse results
python py/parse_hotnet_results.py --time_series "$1" --results_dir "$results"


end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Execution completed in $elapsed_time seconds."
