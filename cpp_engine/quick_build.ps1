#!/usr/bin/env pwsh
$ErrorActionPreference = "Stop"

Write-Host "=== QUICK BUILD ===" -ForegroundColor Green
$SRC = "d:\Chemsi\cpp_engine"
$BUILD = "d:\Chemsi\build-mingw64"

Write-Host "Source: $SRC"
Write-Host "Build: $BUILD"

# Create build dir if needed
if (-not (Test-Path $BUILD)) {
    New-Item -ItemType Directory -Force -Path $BUILD | Out-Null
}

# CMake configure
Set-Location $BUILD
Write-Host "Configuring CMake..." -ForegroundColor Cyan
& cmake.exe -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DEnable_VIS=ON "$SRC"

if ($LASTEXITCODE -ne 0) {
    Write-Host "CMAKE FAILED" -ForegroundColor Red
    exit 1
}

# Build
Write-Host "Building with Ninja..." -ForegroundColor Cyan
& ninja.exe

if ($LASTEXITCODE -ne 0) {
    Write-Host "BUILD FAILED" -ForegroundColor Red
    exit 1
}

Write-Host "=== BUILD SUCCESSFUL ===" -ForegroundColor Green

# Launch
$exe = Join-Path $BUILD "VFEP.exe"
if (Test-Path $exe) {
    Write-Host "Launching visualizer..." -ForegroundColor Green
    & $exe
} else {
    Write-Host "ERROR: VFEP.exe not found at $exe" -ForegroundColor Red
}
