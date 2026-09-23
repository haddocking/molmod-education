#!/bin/bash

file="$1"

# Extract native score from header
native_score=$(awk '/^# native score/ { print $5 }' "$file")

# Extract all unique chains (column 1, skipping header lines)
chains=$(awk 'NR > 10 { print $1 }' "$file" | sort -u)

for chain in $chains; do
    echo "==============="
    echo "== chain $chain"
    
    # Process each chain
    awk -v chain="$chain" -v native="$native_score" '
    $1 == chain {
        score = $5
        if (min == "" || score < min) {
            min = score
            min_mutated = $4
        }
        if (max == "" || score > max) {
            max = score
            max_mutated = $4
        }
    }
    END {
        printf "Mutation of %s leads to the lowest delta haddock score of %.2f\n", min_mutated, min
        printf "Mutation of %s leads to the highest delta haddock score of %.2f\n", max_mutated, max
    }' "$file"
done
echo "==============="
echo "Negative delta haddock score ~ WT presumably stronger binding"
echo "Positive delta haddock score ~ Mutant presumably stronger binding"
