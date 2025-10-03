#!/usr/bin/env bash
set -euo pipefail

INSTALL_DIR=${1:-"$PWD/MAGScoT"}
CONFIG_FILE=${2:-"sunbeam_config.yml"}

echo "Installing MAGScoT into: $INSTALL_DIR"

# Clone if not already present
if [ ! -d "$INSTALL_DIR/.git" ]; then
    git clone https://github.com/ikmb/MAGScoT.git "$INSTALL_DIR"
else
    echo "MAGScoT already present in $INSTALL_DIR. Pulling latest changes..."
    git -C "$INSTALL_DIR" pull
fi

# Update config
if [ -f "$CONFIG_FILE" ]; then
    echo "Updating $CONFIG_FILE with MAGScoT_fp..."
    # remove existing line if present
    grep -v "MAGScoT_fp:" "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" || true
    mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"
    # append the new entry under magscot section
    {
        echo ""
        echo "magscot:"
        echo "  MAGScoT_fp: \"$INSTALL_DIR\""
    } >> "$CONFIG_FILE"
else
    echo "WARNING: $CONFIG_FILE not found. Please add this manually:"
    echo "magscot:"
    echo "  MAGScoT_fp: \"$INSTALL_DIR\""
fi

echo "Done. You can now reference Cfg[\"magscot\"][\"MAGScoT_fp\"] in your rules."
