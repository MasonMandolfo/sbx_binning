#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Usage: vamb_clusters_to_bins.py clusters.tsv > scaffolds2bin.tsv\n")
    sys.exit(1)

clusters = sys.argv[1]

with open(clusters) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        cluster, contig = parts[0], parts[1]
        print(f"{contig}\t{cluster}")
