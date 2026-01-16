@echo off
setlocal

REM Change this only if your MSYS2 is installed elsewhere
set "MSYS=C:\msys64"

if not exist "%MSYS%\usr\bin\bash.exe" (
  echo [ERROR] Cannot find MSYS2 bash at:
  echo   "%MSYS%\usr\bin\bash.exe"
  echo Fix: edit MSYS= in this .bat to your MSYS2 install folder.
  pause
  exit /b 1
)

REM Run the stable bash launcher
"%MSYS%\usr\bin\bash.exe" -lc "bash /d/Chemsi/launch_all.sh"
pause
