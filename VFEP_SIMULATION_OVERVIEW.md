# VFEP Simulation — External Overview

## Status (Jan 23, 2026)
Build pipeline is in place; simulation server and visualization are ready to run once dependencies are resolved.

## What This Is
VFEP is a deterministic fire-suppression simulation and visualization system. It couples a C++ physics/chemistry core with a real-time UI that streams telemetry via gRPC for live inspection and control.

## Quick Start (Runtime)
1. Start VFEP_Sim.exe (gRPC server).
2. Start VFEP.exe (visualization client).

Full run steps: [QUICKSTART.md](QUICKSTART.md).

## Build Workflow (High Level)
1. Tool check for CMake and MinGW64 toolchain.
2. Stop any running VFEP processes to unlock build directory.
3. Clean build directory (or at minimum remove CMake cache).
4. Configure with CHEMSI_BUILD_VIS=ON.
5. Build release binaries.
6. Run NumericIntegrity.exe gate.
7. Launch VFEP.exe (UI).

Build scripts: [launch_all.sh](launch_all.sh) and [cpp_engine/launch_all.bat](cpp_engine/launch_all.bat).

## Core Components
- **Simulation server:** VFEP_Sim.exe — runs the physics/chemistry loop and exposes gRPC telemetry.
- **Visualization client:** VFEP.exe — ImGui/ImPlot UI with 3D scene and live plots.
- **Build scripts:** launch_all.sh and launch_all.bat — standard build/configure/run entry points.

## Key Features (Current Stage)
- Deterministic stepping with telemetry cadence for reproducibility.
- Nozzle automation tied to suppression state.
- Four-tab UI (EXEC, NOZZLE, VIZ, PLOTS) for organized workflows.
- Real-time plots: Temperature, HRR, Exposure, Knockdown, O2.

## Outputs & Artifacts
- Telemetry CSV outputs (e.g., high_fidelity_ml.csv) for analysis.
- Build artifacts in build-mingw64 or build-clean directories.

## Dependencies (Summary)
- CMake, MinGW64 (GCC), Ninja or MinGW Makefiles
- OpenGL/GLFW
- gRPC/Protobuf
- ImGui/ImPlot

## Source Map
- Core simulation: [cpp_engine/src](cpp_engine/src)
- Public API: [cpp_engine/include](cpp_engine/include)
- Visualization: [cpp_engine/vis/main_vis.cpp](cpp_engine/vis/main_vis.cpp)
- Proto definitions: [proto](proto)

## Notes
This document is the single external-facing summary for the current stage. Historical logs and internal verification notes have been archived to reduce noise.
