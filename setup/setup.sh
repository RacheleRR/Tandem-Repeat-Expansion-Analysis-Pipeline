#!/usr/bin/env bash
set -euo pipefail

echo "Installing external tools..."
bash setup/install_tools.sh

echo "Downloading resources..."
bash setup/download_resources.sh

echo "Downloading helper scripts..."
bash setup/download_helpers.sh

echo "Setup complete."
