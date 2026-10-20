#!/bin/bash

set -e

start=$(date +%s.%N)

wolframscript -file example.wl
python3 example.py

end=$(date +%s.%N)

elapsed=$(echo "$end - $start" | bc)
echo "Total execution time: ${elapsed} seconds"
