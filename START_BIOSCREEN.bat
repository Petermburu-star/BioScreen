@echo off
title BioScreen
cd /d "%~dp0"
echo Starting BioScreen - offline biosecurity screening...
python "%~dp0scripts\run_product.py"
pause
