#!/usr/bin/env bash
set -euo pipefail

# Location for MAGScoT clone (inside the extension repo)
MAGSCOt_DIR="$(dirname $0)/../MAGScoT"

# Clone once if missing
if [ ! -d "$MAGSCOt_DIR" ]; then
    git clone --depth 1 https://github.com/ikmb/MAGScoT "$MAGSCOt_DIR"
fi

# Forward all arguments directly to the Rscript
Rscript "$MAGSCOt_DIR/MAGScoT.R" "$@"
