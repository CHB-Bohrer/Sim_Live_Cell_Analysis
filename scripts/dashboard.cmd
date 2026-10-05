@echo off
rem Opens the dashboard in your browser (http://localhost:8501, this PC only). Close this window / Ctrl-C to stop.
call "%~dp0run.cmd" python -m streamlit run "%~dp0..\app\dashboard.py" --server.port 8501 --server.address localhost
