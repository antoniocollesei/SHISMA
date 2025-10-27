import os
import re
import subprocess
import time
import csv
from pathlib import Path

PROJECT_DIR = Path.home() / "Projects" / "SHISMA"
DATA_DIR = PROJECT_DIR / "tests" / "time_test" / "datasets"
WRAPPER = PROJECT_DIR / "SHISMA_main" / "wrapper.sh"
CSV_FILE = PROJECT_DIR / "tests" / "time_test" / "timing_results.csv"

# cd to wrapper directory
os.chdir(WRAPPER.parent)

# Prepare CSV file and write header if not present
write_header = not CSV_FILE.exists()
with open(CSV_FILE, "a", newline="") as csvfile:
  writer = csv.writer(csvfile)
  if write_header:
    writer.writerow(["dataset", "ppi_type", "rows", "cols", "runtime_sec", "timestamp"])

  for datafile in DATA_DIR.glob("subdataset*.csv"):
    prefix = datafile.stem
    ppi = DATA_DIR / f"ppi_{'_'.join(prefix.split('_')[1:-2])}_{prefix.split('_')[-2]}_{prefix.split('_')[-1]}.tsv"

    # Extract rows and cols from filename
    m = re.match(r".*_(full|half|quarter)_(\d+)rows_(\d+)cols", prefix)
    ppi_type = m.group(1) if m else None
    rows, cols = (int(m.group(2)), int(m.group(3))) if m else (None, None)

    print(f"Running {prefix} ...")

    start = time.time()
    subprocess.run([
      str(WRAPPER),
      "--data", str(datafile),
      "--ppi", str(ppi),
      "--ct", "Fakecells",
      "--out", prefix,
      "--nperm", "100",
      "--cores", "25"
    ], check=True)
    runtime = time.time() - start

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    writer.writerow([prefix, ppi_type, rows, cols, runtime, timestamp])

    print(f"{prefix} completed in {runtime:.2f}s")
