# Phase 8 Session Initialization

**Date**: February 2, 2026  
**Phase**: Phase 8 - Advanced Scenarios & Multi-Zone Physics  
**Status**: Ready to start

---

## Session Start Checklist

### ✅ Prerequisites Verified
- [x] Phase 7 Week 1 complete (SensitivityAnalysis, UQ modules)
- [x] 57/57 numeric integrity tests passing
- [x] 4/4 validation scenarios passing
- [x] Clean build (zero warnings/errors)
- [x] PHASE8_STARTUP.md reviewed

### 📋 Phase 8 Objectives
1. **New Scenarios**: Ship Fire, Tunnel Fire, Industrial Compartment (3 scenarios)
2. **Three-Zone Model**: Coupled hot/middle/ambient layers with mass/energy exchange
3. **CFD Interface**: Mock coupling for velocity field import/export
4. **Performance**: Maintain <2s per 60s simulation (currently ~0.3s)
5. **Documentation**: Complete API reference, user guide, technical report

---

## Current System State

### Build Status
```bash
Build: ✅ Clean (MinGW64)
Tests: 57/57 passing
Validation: 4/4 scenarios passing
Warnings: 0
Errors: 0
```

### Validation Results
| Scenario | Error | Status |
|----------|-------|--------|
| ISO 9705 | 4.11% | ✅ PASS |
| NIST | 4.85% | ✅ PASS |
| Suppression | 13.95% | ✅ PASS |
| Stratification | 9.26% | ✅ PASS |

### Phase 7 Modules (Complete)
- ✅ SensitivityAnalysis (parameter sweeps)
- ✅ UncertaintyQuantification (Monte Carlo + LHS)
- ✅ SweepTool (CLI utility)
- ✅ 6 Phase 7 tests (7A1-7A3, 7B1-7B3)

---

## Phase 8 Work Plan (4 Weeks)

### Week 1-2: New Fire Scenarios
**Goal**: Add 3 new validated scenarios

**Week 1:**
- [ ] Day 1-2: Ship Fire scenario (confined, high heat loss)
- [ ] Day 3-4: Tunnel Fire scenario (flow-driven, stratified)
- [ ] Day 5: Testing and validation

**Week 2:**
- [ ] Day 1-2: Industrial Fire scenario (large volume, mixed materials)
- [ ] Day 3-4: Integration with ValidationSuite
- [ ] Day 5: Verify 7/7 scenarios passing

**Deliverables:**
- Ship Fire: <15% error vs IMO SOLAS
- Tunnel Fire: <20% error vs EUREKA 499
- Industrial Fire: <18% error vs ISO 9414

---

### Week 2-3: Three-Zone Model
**Goal**: Implement coupled three-zone physics

**Implementation:**
- [ ] Zone structure (upper/middle/lower with heights)
- [ ] Mass exchange (buoyancy-driven flow)
- [ ] Heat transfer (zone-to-zone + walls)
- [ ] Species transport (advection + diffusion)
- [ ] Validation against ISO 9705

**Deliverables:**
- ThreeZoneFireModel operational
- ISO 9705 with three-zone: 1000-1050K (vs 981K single-zone)
- 3 numeric tests (8A1-8A3)

---

### Week 3-4: CFD Interface & Documentation
**Goal**: Define CFD coupling and complete documentation

**CFD Interface:**
- [ ] VTK import/export
- [ ] Trilinear interpolation
- [ ] Mock CFD generator
- [ ] Comparison workflow
- [ ] 2 numeric tests (8B1-8B2)

**Documentation:**
- [ ] Phase8_API_Reference.md
- [ ] Phase8_User_Guide.md
- [ ] Phase8_Technical_Report.md
- [ ] Phase8_Scenarios_Catalog.md

**Deliverables:**
- CFDCoupler interface complete
- 4 documentation files published
- 62/62 tests passing
- 7/7 scenarios passing

---

## Implementation Order

### Priority 1 (Week 1): Ship Fire
```cpp
// ValidationSuite.cpp additions
void validateShipFire() {
  Simulation sim;
  sim.resetToScenario({
    .volume_m3 = 100.0,
    .area_m2 = 150.0,
    .h_W_m2K = 15.0,  // Aluminum
    .ACH = 0.8,
    .pyrolysis_kgps = 0.10,
    .heatRelease_J_per_molFuel = 150e3
  });
  
  // Run and validate
  // Target: 800-1000K peak temperature
}
```

### Priority 2 (Week 1): Tunnel Fire
```cpp
void validateTunnelFire() {
  Simulation sim;
  sim.resetToScenario({
    .volume_m3 = 5000.0,
    .area_m2 = 1000.0,
    .h_W_m2K = 8.0,   // Concrete
    .ACH = 5.0,
    .pyrolysis_kgps = 0.20,
    .heatRelease_J_per_molFuel = 450e3
  });
  
  // Target: 500-2000 kW HRR plateau
}
```

### Priority 3 (Week 2): Industrial Fire
```cpp
void validateIndustrialFire() {
  Simulation sim;
  sim.resetToScenario({
    .volume_m3 = 2500.0,
    .area_m2 = 800.0,
    .h_W_m2K = 5.0,   // Steel
    .ACH = 3.0,
    .pyrolysis_kgps = 0.35,
    .heatRelease_J_per_molFuel = 200e3
  });
  
  // Target: >500K ceiling temperature
}
```

### Priority 4 (Week 2-3): Three-Zone Model
```cpp
// New file: ThreeZoneModel.h
class ThreeZoneFireModel {
  Zone upper_, middle_, lower_;
  
  void step(double dt, double HRR_W, double cooling_W, double ACH);
  void updateZoneBoundaries();
  void updateMassExchange();
  void updateHeatTransfer();
};
```

### Priority 5 (Week 3): CFD Interface
```cpp
// New file: CFDInterface.h
class CFDCoupler {
  void importVelocityField(const std::string& vtk_file);
  void exportResults(const std::string& output_vtk);
  ComparisonStats compareTemperature(const std::string& cfd_file);
};
```

---

## Testing Strategy

### Numeric Integrity Tests (Phase 8)

**Total Target**: 62 tests (57 existing + 5 Phase 8)

**New Phase 8 Tests:**
- 8A1: Three-Zone Basic Initialization
- 8A2: Three-Zone Mass Conservation
- 8A3: Three-Zone Energy Balance
- 8B1: CFD Import Basic
- 8B2: CFD Export and Compare

### Validation Scenarios

**Total Target**: 7 scenarios

**Existing (4):**
- ✅ ISO 9705 (4.11%)
- ✅ NIST (4.85%)
- ✅ Suppression (13.95%)
- ✅ Stratification (9.26%)

**New (3):**
- 🆕 Ship Fire (<15%)
- 🆕 Tunnel Fire (<20%)
- 🆕 Industrial Fire (<18%)

---

## File Structure (Phase 8 Additions)

```
cpp_engine/
├── include/
│   ├── ThreeZoneModel.h        (NEW)
│   └── CFDInterface.h          (NEW)
├── src/
│   ├── ThreeZoneModel.cpp      (NEW)
│   └── CFDInterface.cpp        (NEW)
├── tools/
│   └── GenerateMockCFD.cpp     (NEW)
└── tests/
    └── TestNumericIntegrity.cpp (UPDATE: add 8A1-8A3, 8B1-8B2)

ValidationSuite.cpp (UPDATE: add 3 scenarios)

Documentation:
├── PHASE8_API_Reference.md     (NEW)
├── PHASE8_User_Guide.md        (NEW)
├── PHASE8_Technical_Report.md  (NEW)
└── PHASE8_Scenarios_Catalog.md (NEW)
```

---

## Daily Workflow

### Morning Routine
1. Verify baseline: `.\NumericIntegrity.exe` (should show 57/57)
2. Verify validation: `.\ValidationSuite.exe` (should show 4/4)
3. Check build: `cmake --build build-mingw64 --config Release`

### Development Cycle
1. Implement feature (scenario, zone model, CFD interface)
2. Add numeric tests immediately
3. Build and verify tests pass
4. Validate against literature
5. Document in appropriate .md file

### End of Day
1. Run full test suite (NumericIntegrity + ValidationSuite)
2. Commit working code
3. Update progress in this document

---

## Success Metrics (Phase 8 Complete)

**Required for Phase 8 sign-off:**

1. ✅ **62/62** numeric integrity tests passing
2. ✅ **7/7** validation scenarios passing
3. ✅ **Zero** compiler warnings/errors
4. ✅ **ThreeZoneFireModel** operational and validated
5. ✅ **CFDCoupler** interface complete
6. ✅ **4 documentation files** published
7. ✅ **Performance** <2s per 60s simulation maintained

---

## Resources

### Primary Documents
- [PHASE8_STARTUP.md](PHASE8_STARTUP.md) - Complete work plan
- [PHASE7_INTEGRITY_REPORT.md](PHASE7_INTEGRITY_REPORT.md) - Baseline state
- [readme.md](readme.md) - Project overview

### Literature References
- **IMO SOLAS** - Ship fire standards
- **EUREKA 499** - Tunnel fire experiments
- **Memorial Tunnel** - Ventilation test program
- **ISO 9414** - Gross calorific potential
- **ASTM E603** - Room fire test

### Build Commands
```bash
# Configure
cmake -S . -B build-mingw64 -G "MinGW Makefiles"

# Build
cmake --build build-mingw64 --config Release

# Test
cd build-mingw64
.\NumericIntegrity.exe
.\ValidationSuite.exe
```

---

## Common Issues & Solutions

### Issue: New scenario fails validation
**Solution**: Check parameter ranges, verify literature source, adjust tolerance if physics-justified

### Issue: Three-zone model unstable
**Solution**: Reduce timestep, add mass/energy guards, verify zone boundary logic

### Issue: CFD import fails
**Solution**: Verify VTK format, check file paths, validate grid dimensions

### Issue: Performance degrades
**Solution**: Profile with `gprof`, identify hotspots, optimize selectively

---

## Communication Protocol

### Daily Updates
Report:
- Tests passing (X/62)
- Scenarios passing (Y/7)
- Blockers (if any)
- Next day plan

### Week Milestones
- Week 1: 3 new scenarios
- Week 2: ValidationSuite 7/7
- Week 3: Three-zone + CFD
- Week 4: Documentation complete

---

## Let's Build Phase 8! 🚀

**First task**: Begin Ship Fire scenario implementation in ValidationSuite.cpp

**Expected time**: 2-3 hours for basic scenario + validation

**Success criteria**: Ship Fire passes with <15% error vs IMO SOLAS data
