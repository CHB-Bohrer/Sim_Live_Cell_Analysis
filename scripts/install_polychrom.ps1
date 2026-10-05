# Install polychrom from GitHub WITHOUT its optional Cython extension (_polymer_math, used only for
# knot simplification). Avoids needing MSVC Build Tools on Windows.
# Usage (inside the activated `simlive` env):  powershell -File scripts/install_polychrom.ps1
$ErrorActionPreference = "Stop"
$tmp = Join-Path $env:TEMP "polychrom_src"
if (Test-Path $tmp) { Remove-Item -Recurse -Force $tmp }
git clone --depth 1 https://github.com/open2c/polychrom.git $tmp
Set-Content -Path (Join-Path $tmp "setup.py") -Value "from setuptools import setup`nsetup()`n"
python -m pip install $tmp
