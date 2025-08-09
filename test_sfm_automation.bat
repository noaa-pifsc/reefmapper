@echo off
REM Set the path to Metashape python interpreter
set METASHAPE_PATH="C:\Program Files\Agisoft\Metashape Pro\python\python.exe"

REM Set the path to your Python script
set SCRIPT_PATH="C:\Users\PICHLMRUser\Desktop\reefmapper\SfMBatchProcess_db.py"

REM Call Metashape with your script
%METASHAPE_PATH% %SCRIPT_PATH%

pause