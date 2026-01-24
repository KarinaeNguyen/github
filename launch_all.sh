#!/usr/bin/env bash
set -euo pipefail

ROOT=/d/Chemsi
SRC="$ROOT/cpp_engine"
BUILD="$ROOT/build-mingw64"

echo "=== TOOL CHECK: CMake ==="

# Prefer MSYS2 cmake if available
if command -v cmake >/dev/null 2>&1; then
  CMAKE="cmake"
else
  # Fallback: try common Windows CMake locations
  if [ -x "/c/Program Files/CMake/bin/cmake.exe" ]; then
    CMAKE="/c/Program Files/CMake/bin/cmake.exe"
  elif [ -x "/c/Program Files (x86)/CMake/bin/cmake.exe" ]; then
    CMAKE="/c/Program Files (x86)/CMake/bin/cmake.exe"
  else
    echo "ERROR: cmake not found in this MINGW64 environment."
    echo ""
    echo "Fix option A (recommended): install CMake in MSYS2:"
    echo "  pacman -S --needed mingw-w64-x86_64-cmake"
    echo ""
    echo "Fix option B: install CMake for Windows and ensure it is in:"
    echo "  C:\\Program Files\\CMake\\bin\\cmake.exe"
    exit 127
  fi
fi

echo "Using CMake: $CMAKE"
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
"$CMAKE" -S "$SRC" -B "$BUILD" -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCHEMSI_BUILD_VIS=ON

echo "=== BUILD ==="
"$CMAKE" --build "$BUILD" -j

echo "=== NUMERIC INTEGRITY GATE ==="
"$BUILD/NumericIntegrity.exe"

echo "=== LAUNCH UI ==="
if [[ ! -f "$BUILD/VFEP.exe" ]]; then
  echo "ERROR: $BUILD/VFEP.exe not found (visualizer not built)."
  echo "Build folder contents:"
  ls -la "$BUILD" || true
  exit 1
fi

# Ensure MinGW runtime DLLs are found
export PATH=/c/msys64/mingw64/bin:$PATH

# Convert VFEP path to Windows path for start
# (removed) launching directly; no Windows path needed

# Launch detached so UI stays open
exec "$BUILD/VFEP.exe" --calib

echo "UI launched."
