@echo off
setlocal

set "MSYS=C:\msys64"

if not exist "%MSYS%\usr\bin\bash.exe" (
  echo [ERROR] Cannot find MSYS2 bash at: "%MSYS%\usr\bin\bash.exe"
  pause
  exit /b 1
)

REM Force MinGW64 environment PATH so cmake/g++/mingw32-make are visible
"%MSYS%\usr\bin\bash.exe" -lc "export MSYSTEM=MINGW64; export PATH=/mingw64/bin:/usr/bin:$PATH; bash /d/Chemsi/launch_all.sh"
pause
