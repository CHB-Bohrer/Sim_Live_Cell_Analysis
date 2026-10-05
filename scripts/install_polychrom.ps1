# Install polychrom from GitHub WITHOUT its optional Cython extension (_polymer_math, used only for
# knot simplification). Avoids needing MSVC Build Tools on Windows.
# Pinned to the exact commit the project was validated with (see environments/PINNED.md), so a rebuild gives the same
# polychrom even if upstream has changed.
# Usage (inside the activated `simlive` env):  powershell -File scripts/install_polychrom.ps1
#   to try a different commit (in a TEST environment first):  ... install_polychrom.ps1 -Commit <sha>
param([string]$Commit = "11a870cac8a3b168a2e0e11995d2e899c3ceb657")
$ErrorActionPreference = "Stop"
$tmp = Join-Path $env:TEMP "polychrom_src"
if (Test-Path $tmp) { Remove-Item -Recurse -Force $tmp }
git clone https://github.com/open2c/polychrom.git $tmp
git -C $tmp checkout $Commit
Set-Content -Path (Join-Path $tmp "setup.py") -Value "from setuptools import setup`nsetup()`n"
python -m pip install $tmp
