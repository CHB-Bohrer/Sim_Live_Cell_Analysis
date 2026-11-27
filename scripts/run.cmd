@echo off
rem Run any command inside the `simlive` conda env, from any prompt (no activation needed).
rem   scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\demo
rem   scripts\run.cmd pytest -m gpu
rem Simulation outputs (runs, sweeps, chromatin library) live in SIMLIVE_DATA, outside the git repo:
rem by default the folder `simulation_data` next to the repo folder (Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data),
rem else %USERPROFILE%\SimLiveData (the old place; on this PC a link to the new one).
if not defined SIMLIVE_DATA (
  if exist "%~dp0..\..\simulation_data\" (
    for %%I in ("%~dp0..\..\simulation_data") do set "SIMLIVE_DATA=%%~fI"
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
