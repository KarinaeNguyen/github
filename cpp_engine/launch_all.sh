#!/usr/bin/env bash
set -euo pipefail

ROOT=/d/Chemsi
SRC="$ROOT/cpp_engine"
BUILD="$ROOT/build-mingw64"

echo "=== CLEAN ==="
rm -rf "$BUILD"

echo "=== CONFIGURE (SRC=$SRC, BUILD=$BUILD, VIS ON) ==="
cmake -S "$SRC" -B "$BUILD" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCHEMSI_BUILD_VIS=ON

echo "=== BUILD ==="
cmake --build "$BUILD" -j

echo "=== NUMERIC INTEGRITY GATE ==="
"$BUILD/NumericIntegrity.exe"

echo "=== LAUNCH UI ==="
if [[ ! -f "$BUILD/VFEP.exe" ]]; then
  echo "ERROR: $BUILD/VFEP.exe not found."
  echo "Contents of build folder:"
  ls -ლა "$BUILD" || ls -la "$BUILD"
  exit 1
fi

# Ensure MinGW runtime DLLs are found
export PATH=/c/msys64/mingw64/bin:$PATH

# Convert to Windows path for cmd.exe / start
VFEP_WIN="$(cygpath -w "$BUILD/VFEP.exe")"

# Launch detached so the UI stays open
cmd.exe /c start "" "%VFEP_WIN%" --calib

echo "UI launched: $BUILD/VFEP.exe"
