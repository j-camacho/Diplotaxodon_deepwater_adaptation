#!/bin/bash

IN_DIR=path_to_dir/HTSeq
OUT_DIR=path_to_dir/HTSeq/counts
MT=path_to_dir/metadata
PATTERN="Aligned.sortedByCoord.out.count.txt"


for file in "$IN_DIR"/*"$PATTERN"; do
    sample_name=$(basename "$file" | sed "s/$PATTERN//")
    echo "$sample_name"
    awk '{print $2}' "$file" > "$OUT_DIR"/"${sample_name}".count
done

awk '{print $1}' "$file" > "$OUT_DIR/geneids.txt"

paste "$OUT_DIR/geneids.txt" "$OUT_DIR"/*.count > "$OUT_DIR/counts_allsamples.tmp.out"

cat <(cat "$MT/samples.txt" | sort | paste -s) "$OUT_DIR/counts_allsamples.tmp.out" > "$OUT_DIR/06_rnaseq_count-matrix_htseq-count.txt"

rm "$OUT_DIR"/*.count "$OUT_DIR/counts_allsamples.tmp.out"
