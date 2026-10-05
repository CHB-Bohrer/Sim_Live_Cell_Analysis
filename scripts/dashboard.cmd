@echo off
rem Starts the dashboard and opens it in your browser (http://127.0.0.1:8501, this PC only).
rem KEEP THIS WINDOW OPEN while you use the dashboard. Close it (or Ctrl-C) to stop.
cd /d "%~dp0.."
echo Starting dashboard on http://127.0.0.1:8501 ... your browser will open in about 10 seconds.
echo KEEP THIS WINDOW OPEN. Close it to stop the dashboard.
echo.
start "" /b cmd /c "ping -n 11 127.0.0.1 >nul & start "" http://127.0.0.1:8501"
call "%~dp0run.cmd" python -m streamlit run app\dashboard.py --server.port 8501 --server.address 127.0.0.1 --server.headless true --browser.gatherUsageStats false
echo.
echo The dashboard has stopped. If this was not intentional, the messages above explain why.
pause
