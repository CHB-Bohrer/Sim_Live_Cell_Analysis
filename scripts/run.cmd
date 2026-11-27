@echo off
rem Run any command inside the `simlive` conda env, from any prompt (no activation needed).
rem   scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\demo
rem   scripts\run.cmd pytest -m gpu
rem Simulation outputs (runs, sweeps, chromatin library) live in SIMLIVE_DATA, outside the git repo:
rem first match wins (src\simlive\io\runs.py DATA_ROOT uses the same order): a folder `simulation_data` next to the
rem repo folder; the owner's Dropbox data folder %USERPROFILE%\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data
rem (the repo itself lives in %USERPROFILE%\Projects\Sim_Live_Cell_Analysis since 2026-11-27); else %USERPROFILE%\SimLiveData
rem (the old place; on the owner's PC a link to the Dropbox folder).
if not defined SIMLIVE_DATA (
  if exist "%~dp0..\..\simulation_data\" (
    for %%I in ("%~dp0..\..\simulation_data") do set "SIMLIVE_DATA=%%~fI"
  ) else if exist "%USERPROFILE%\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data\" (
    set "SIMLIVE_DATA=%USERPROFILE%\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data"
  ) else (
    set "SIMLIVE_DATA=%USERPROFILE%\SimLiveData"
  )
)
if not exist "%USERPROFILE%\miniforge3\condabin\conda.bat" (
  echo Miniforge not found at %USERPROFILE%\miniforge3. See README.md "Setup".
  exit /b 1
)
call "%USERPROFILE%\miniforge3\condabin\conda.bat" activate simlive || exit /b 1
%*
