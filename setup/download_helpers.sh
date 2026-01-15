#!/usr/bin/env bash
set -euo pipefail

mkdir -p helper
cd helper

if [ ! -d "Fazal2020Scripts" ]; then
    git clone https://github.com/ZuchnerLab/Fazal2020Scripts.git
    cd Fazal2020Scripts
    chmod +x *
fi
cd ..
if [ ! -d "BTlib" ]; then
git clone https://github.com/bjtrost/tandem-repeat-expansions-in-ASD.git
cd tandem-repeat-expansions-in-ASD
chmod +x *
cp -r BTlib ..
cp BTlib.py ..
cd ..
rm -rf tandem-repeat-expansions-in-ASD
fi 

patch -N < BTlib_docx_imports.patch || true
