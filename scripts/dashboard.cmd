@echo off
rem Starts the dashboard and opens it in your browser (http://localhost:8501, this PC only).
rem Leave this window open while you use it; close it or press Ctrl-C to stop.
echo Starting dashboard... your browser will open in a few seconds. Keep this window open.
start "" /b cmd /c "timeout /t 8 /nobreak >nul & start "" http://localhost:8501"
call "%~dp0run.cmd" python -m streamlit run "%~dp0..\app\dashboard.py" --server.port 8501 --server.address localhost --server.headless true --browser.gatherUsageStats false
