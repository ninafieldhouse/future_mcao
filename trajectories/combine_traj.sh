#!/bin/bash

# Usage: ./combine_traj.sh output.txt input1.txt input2.txt input3.txt ...

output="$1"
shift

first=1
for file in "$@"; do
    if [ $first -eq 1 ]; then
        # Copy everything (including header) from the first file
        cat "$file" > "$output"
        first=0
    else
        # Skip the first 4 lines (header) for subsequent files
        tail -n +5 "$file" >> "$output"
    fi
done

echo "Combined files into $output"