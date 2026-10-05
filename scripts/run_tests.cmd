@echo off
rem Runs ALL tests in two separate processes (GPU tests and CPU tests cannot share one process on Windows).
call "%~dp0run.cmd" pytest -m gpu -q || exit /b 1
call "%~dp0run.cmd" pytest -q
