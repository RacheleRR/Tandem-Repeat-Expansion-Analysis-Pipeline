#!/usr/bin/env bash
set -euo pipefail

mkdir -p program
cd program

if [ ! -d "ehdn" ]; then
    echo "Installing ExpansionHunterDenovo v0.9.0"

    wget https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

    tar -xzf ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

    mv ExpansionHunterDenovo-v0.9.0-linux_x86_64 ehdn

    rm ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz
fi

if [ ! -d "liftOver" ]; then
     wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver   liftOver 
        chmod +x liftOver
fi

if [ ! -d "rex" ]; then
    git clone https://github.com/ZuchnerLab/RExPRT.git RExPRT
    cd RExPRT/RExPRT
    chmod +x RExPRT.sh
fi
