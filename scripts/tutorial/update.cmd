@echo off
rem Keeps the narrated tutorial in step with the dashboard. Details: docs\tutorial\UPDATING.md
rem   update.cmd            what is out of date + a check of every sentence
rem   update.cmd all        re-take stale screenshots, record edited sentences, check
setlocal
set "PY=%TUTORIAL_PY%"
if "%PY%"=="" set "PY=%USERPROFILE%\shotenv_tmp\v\Scripts\python.exe"
if not exist "%PY%" (
  echo No Python with playwright at %PY%
  echo Make one:  python -m venv %USERPROFILE%\shotenv_tmp\v   then   %USERPROFILE%\shotenv_tmp\v\Scripts\pip install playwright imageio-ffmpeg
  echo or set TUTORIAL_PY to a python.exe that has playwright.
  exit /b 1
)
"%PY%" "%~dp0update_tutorial.py" %*
