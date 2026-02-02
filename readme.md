# CHEMSI / VFEP Simulation

**Current Project Phase**: Phase 7 - ✅ **IN PROGRESS (Sensitivity Analysis & UQ)**

## Quick Navigation

### 📋 Start Here
- **New to project?** → [VFEP_SIMULATION_OVERVIEW.md](VFEP_SIMULATION_OVERVIEW.md)
- **Run quickly?** → [QUICKSTART.md](QUICKSTART.md)
- **Phase 7 status?** → [PHASE7_INTEGRITY_REPORT.md](PHASE7_INTEGRITY_REPORT.md)

### 🎯 Current Status
- **Validation Rate**: 100% (4/4 scenarios passing)
- **Numeric Tests**: 57/57 passing (including Phase 7 tests)
- **Build Status**: ✅ Clean, no warnings
- **Phase 7 Modules**: SensitivityAnalysis, UncertaintyQuantification

### 📊 Phase 7 Progress
- ✅ SensitivityAnalysis module (parameter sweeps)
- ✅ UncertaintyQuantification module (Monte Carlo + LHS)
- ✅ 6 numeric integrity tests (7A1-7A3, 7B1-7B3)
- ✅ SweepTool CLI utility
- ⏳ New scenarios (Ship Fire, Tunnel Fire, Industrial)
- ⏳ Three-zone model implementation
- ⏳ CFD coupling interface (mock)

### 📊 Validation Results
| Scenario | Predicted | Target | Error | Status |
|----------|-----------|--------|-------|--------|
| ISO 9705 | 981K | 1023K ±50K | 4.11% | ✅ PASS |
| NIST Data Center | 71.4 kW | 75 kW | 4.85% | ✅ PASS |
| Suppression | 79.77% | 60-80% | 13.95% | ✅ PASS |
| Stratification | 272K ΔT | 200-400K | 9.26% | ✅ PASS |

## Key Paths
- Core simulation: [cpp_engine/src](cpp_engine/src)
- Public API: [cpp_engine/include](cpp_engine/include)
- Visualization: [cpp_engine/vis/main_vis.cpp](cpp_engine/vis/main_vis.cpp)
- Build scripts: [launch_all.sh](launch_all.sh) and [cpp_engine/launch_all.bat](cpp_engine/launch_all.bat)
- Protos: [proto](proto)

## Phase History
- **Phase 1-2**: Foundation (combustion, thermodynamics)
- **Phase 3-4**: Integration (suppression, ventilation)
- **Phase 5**: Initial calibration (NIST baseline)
- **Phase 6**: Multi-scenario validation (✅ COMPLETE)
- **Phase 7**: Advanced validation & UQ (✅ WEEK 1 COMPLETE)

Historical logs and verification notes are archived to keep documentation concise.
