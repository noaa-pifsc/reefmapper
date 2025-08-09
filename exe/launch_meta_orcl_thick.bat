@echo off
REM === Set Oracle Instant Client path ===
set ORACLE_CLIENT_PATH=C:\oracle\instantclient_21_19
set PATH=%ORACLE_CLIENT_PATH%;%PATH%
start "" "C:\Program Files\Agisoft\Metashape Pro\metashape.exe"

pause
