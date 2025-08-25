cd ~/Projects/SHISMA/SHISMA_main/

for datafile in ~/Projects/SHISMA/tests/random_data_test/datasets/samples*.csv; do
  base=$(basename "$datafile")
  prefix="${base%.csv}" 
  ppi="~/Projects/SHISMA/tests/random_data_test/datasets/${prefix}_ppi.tsv"
  ./wrapper.sh --data "$datafile" --ppi "$ppi" --ct Bcell --out "$prefix" --nperm 100 --cores 25
done