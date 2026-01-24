#!/usr/bin/env bash
set -euo pipefail

ROOT=/d/Chemsi
SRC="$ROOT/cpp_engine"
BUILD="$ROOT/build-mingw64"

echo "=== TOOL CHECK (must be MinGW64) ==="
echo "MSYSTEM=${MSYSTEM:-"(unset)"}"
command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found on PATH"; exit 127; }
command -v g++  >/dev/null 2>&1 || { echo "ERROR: g++ not found on PATH"; exit 127; }
command -v mingw32-make >/dev/null 2>&1 || { echo "ERROR: mingw32-make not found on PATH"; exit 127; }

echo "cmake:        $(command -v cmake)"
echo "g++:          $(command -v g++)"
echo "mingw32-make: $(command -v mingw32-make)"

echo "SRC:   $SRC"
echo "BUILD: $BUILD"

echo "=== STOP RUNNING PROCESSES ==="
# Gracefully stop any running VFEP processes to unlock build directory
taskkill //F //IM VFEP_Sim.exe 2>/dev/null || echo "  (no VFEP_Sim.exe running)"
taskkill //F //IM VFEP.exe 2>/dev/null || echo "  (no VFEP.exe running)"
sleep 1

echo "=== CLEAN ==="
# Try to clean, but continue if directory is busy
if rm -rf "$BUILD" 2>/dev/null; then
  echo "  Build directory cleaned"
else
  echo "  WARNING: Could not fully clean $BUILD (may be in use)"
  echo "  Removing only CMakeCache.txt to force reconfigure..."
  rm -f "$BUILD/CMakeCache.txt" || true
fi

echo "=== CONFIGURE (VIS ON) ==="
cmake -S "$SRC" -B "$BUILD" -G "MinGW Makefiles" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCHEMSI_BUILD_VIS=ON ^
  -DCMAKE_MAKE_PROGRAM=mingw32-make

echo "=== BUILD ==="
cmake --build "$BUILD" -j

echo "=== NUMERIC INTEGRITY GATE ==="
"$BUILD/NumericIntegrity.exe"

echo "=== LAUNCH UI ==="
if [[ ! -f "$BUILD/VFEP.exe" ]]; then
  echo "ERROR: $BUILD/VFEP.exe not found (visualizer not built)."
  ls -la "$BUILD" || true
  exit 1
fi

# Ensure MinGW runtime DLLs are found
export PATH=/mingw64/bin:$PATH

# Run directly (most reliable). This will open the UI window.
exec "$BUILD/VFEP.exe" --calib

# Ensure MinGW runtime DLLs are found
export PATH=/mingw64/bin:$PATH

VFEP_WIN="$(cygpath -w "$BUILD/VFEP.exe")"
cmd.exe /c start "" "%VFEP_WIN%" --calib

echo "UI launched."

exec "$BUILD/VFEP.exe" --calib