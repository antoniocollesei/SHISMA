cd ~/Projects/SHISMA/SHISMA_main/

for datafile in ~/Projects/SHISMA/tests/injection_test/datasets/synthetic*.csv; do
  base=$(basename "$datafile" .csv)
  ppi="~/Projects/SHISMA/tests/injection_test/datasets/${base}_ppi.tsv"
  out="${base}"
  ./wrapper_from_borf.sh --data "$datafile" --ppi "$ppi" --ct Fakecells --out "$out" --nperm 10 --cores 25
done