#!/usr/bin/env bash
# Install Python dependencies
pip install -r requirements.txt

# Download AutoDock Vina for Linux
echo "Downloading AutoDock Vina (Linux) for Render environment..."
curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -o scripts/vina_linux
chmod +x scripts/vina_linux
