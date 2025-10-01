#!/usr/bin/env bash
# Simple wrapper so Snakemake can just call `magscot`
# Uses the installed MAGScoT.R inside the pip package.

MAGS_DIR=$(python -m pip show MAGScoT | grep Location | awk '{print $2}')/MAGScoT
Rscript $MAGS_DIR/MAGScoT.R "$@"
