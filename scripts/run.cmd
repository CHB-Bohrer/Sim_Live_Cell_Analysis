@echo off
rem Run any command inside the `simlive` conda env, from any prompt (no activation needed).
rem   scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\demo
rem   scripts\run.cmd pytest -m gpu
rem Simulation outputs (runs, sweeps, chromatin library) live in SIMLIVE_DATA, outside the OneDrive-synced repo.
if not defined SIMLIVE_DATA set "SIMLIVE_DATA=%USERPROFILE%\SimLiveData"
if not exist "%USERPROFILE%\miniforge3\condabin\conda.bat" (
  echo Miniforge not found at %USERPROFILE%\miniforge3. See README.md "Setup".
  exit /b 1
)
call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate simlive || exit /b 1
%*
