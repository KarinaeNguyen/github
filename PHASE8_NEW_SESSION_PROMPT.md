# VFEP Phase 8 Development Session - Initialization Prompt

**Copy this entire document and paste into a new conversation to initialize Phase 8 work.**

---

## Project Context

I'm working on VFEP (Variable Fire Environment Physics), a deterministic fire-suppression simulation system written in C++ with a real-time visualization UI. The project is currently at **Phase 8: Advanced Scenarios & Multi-Zone Physics**.

### Current Project State (February 2, 2026)

**Build Status:**
- Platform: Windows with MinGW64 toolchain
- CMake build system (CMakeLists.txt in cpp_engine/)
- Build directory: `d:\Chemsi\build-mingw64`
- All targets compile cleanly (zero warnings, zero errors)

**Test Status:**
- ✅ **57/57 numeric integrity tests** passing (`NumericIntegrity.exe`)
- ✅ **4/4 validation scenarios** passing (`ValidationSuite.exe`)
- All scenarios within ±15% literature bounds

**Validation Results:**
| Scenario | Predicted | Literature | Error | Status |
|----------|-----------|------------|-------|--------|
| ISO 9705 Room Corner | 981K | 973-1073K | 4.11% | ✅ PASS |
| NIST Data Center | 71.4 kW | 60-90 kW | 4.85% | ✅ PASS |
| Suppression Effectiveness | 79.77% | 60-80% | 13.95% | ✅ PASS |
| Temperature Stratification | 272K ΔT | 200-400K | 9.26% | ✅ PASS |

---

## Phase History & Accomplishments

**Phase 1-6:** Core physics implementation and validation
- Single-zone fire model with combustion chemistry
- Suppression system (inhibitor + inert gas)
- Ventilation (ACH-based mole/energy exchange)
- Reactor thermodynamics with Cp-weighted energy balance
- 4 validated scenarios against NIST/SFPE literature

**Phase 7 (Week 1 Complete):**
- ✅ **SensitivityAnalysis module** - Parameter sweeps (heat release, wall loss, geometry, pyrolysis)
- ✅ **UncertaintyQuantification module** - Monte Carlo with Latin Hypercube Sampling
- ✅ **SweepTool CLI utility** - Batch parameter studies
- ✅ **6 new numeric tests** (7A1-7A3 for sensitivity, 7B1-7B3 for UQ)
- ✅ Updated to 57/57 tests total

**Phase 8 (Current - Week 1 Starting):**
Ready to implement advanced scenarios and multi-zone physics.

---

## Phase 8 Objectives (4 Weeks)

### Week 1-2: New Fire Scenarios (Priority)
Add 3 new validated scenarios to expand beyond data center fires:

1. **Ship Fire** (Confined Compartment)
   - Geometry: 100 m³ volume, 150 m² area
   - High heat loss: h_W = 15 W/m²K (aluminum structure)
   - Limited ventilation: ACH = 0.8
   - Cellulose fuel: 150 kJ/mol heat release, 0.10 kg/s pyrolysis
   - Target: 800-1000K peak temperature
   - Literature: IMO SOLAS fire testing standards
   - Validation tolerance: <15% error

2. **Tunnel Fire** (Flow-Driven)
   - Geometry: 5000 m³ volume, 1000 m² area
   - Moderate heat loss: h_W = 8 W/m²K (concrete)
   - Strong ventilation: ACH = 5.0 (equivalent to 1 m/s flow)
   - Hydrocarbon fuel: 450 kJ/mol, 0.20 kg/s pyrolysis
   - Target: 500-2000 kW HRR plateau
   - Literature: EUREKA 499, Memorial Tunnel experiments
   - Validation tolerance: <20% error

3. **Industrial Compartment** (Warehouse)
   - Geometry: 2500 m³ volume, 800 m² area
   - Low heat loss: h_W = 5 W/m²K (steel radiates)
   - Moderate ventilation: ACH = 3.0
   - Mixed fuel: 200 kJ/mol, 0.35 kg/s pyrolysis
   - Target: >500K ceiling temperature (sprinkler activation)
   - Literature: ISO 9414, ASTM E603
   - Validation tolerance: <18% error

### Week 2-3: Three-Zone Model
Implement true three-zone fire model (upper hot layer, middle transition, lower ambient) with:
- Mass exchange (buoyancy-driven inter-zone flow)
- Heat transfer (zone-to-zone + wall losses)
- Species transport (advection + diffusion)
- Validation against ISO 9705 (expect 1000-1050K vs 981K current)

### Week 3-4: CFD Interface & Documentation
- Mock CFD coupling (VTK import/export)
- Complete documentation suite (4 files)
- Performance verification (<2s per 60s simulation - already at ~0.3s)

**Target Completion:** 62/62 tests, 7/7 scenarios, ThreeZoneModel operational, CFDCoupler defined

---

## Codebase Architecture

### Core Simulation Files (Phase 6 - LOCKED, Do Not Modify)
```
cpp_engine/src/
├── Simulation.cpp        - Main simulation loop, telemetry, command handling
├── Reactor.cpp           - Thermodynamics, heat capacity, temperature updates
├── chemistry.cpp         - Combustion kinetics (Arrhenius + pilot floor)
├── Suppression.cpp       - VFEP actuator, inhibitor/inert injection
├── Ventilation.cpp       - ACH-based mole exchange, Cp-weighted energy
├── LiIonRunaway.cpp      - Thermal runaway heat source
└── ValidationSuite.cpp   - Literature-based scenario validation

cpp_engine/include/
├── Simulation.h
├── Reactor.h
├── Chemistry.h
├── Suppression.h
├── Ventilation.h
└── Constants.h
```

### Phase 7 Modules (Complete, Stable)
```
cpp_engine/src/
├── SensitivityAnalysis.cpp    - Parameter sweep framework
└── UncertaintyQuantification.cpp - Monte Carlo + LHS

cpp_engine/include/
├── SensitivityAnalysis.h
└── UncertaintyQuantification.h

cpp_engine/tools/
└── SweepTool.cpp - CLI for batch parameter studies

cpp_engine/tests/
└── TestNumericIntegrity.cpp - 57 tests (includes 7A1-7A3, 7B1-7B3)
```

### Phase 8 Files (To Be Created)
```
cpp_engine/include/
├── ThreeZoneModel.h      - NEW: Three-zone fire model interface
└── CFDInterface.h        - NEW: VTK import/export for CFD coupling

cpp_engine/src/
├── ThreeZoneModel.cpp    - NEW: Zone coupling implementation
└── CFDInterface.cpp      - NEW: File I/O and interpolation

cpp_engine/tools/
└── GenerateMockCFD.cpp   - NEW: Synthetic CFD data generator

Documentation (NEW):
├── PHASE8_API_Reference.md
├── PHASE8_User_Guide.md
├── PHASE8_Technical_Report.md
└── PHASE8_Scenarios_Catalog.md
```

---

## Key Technical Details

### Simulation API (resetToScenario pattern)
```cpp
// Example: ISO 9705 scenario
sim.resetToDataCenterRackScenario(); // Built-in

// Custom scenario (what we'll use for Phase 8)
Simulation sim;
sim.resetToScenario({
  .volume_m3 = 100.0,
  .area_m2 = 150.0,
  .h_W_m2K = 15.0,
  .ACH = 0.8,
  .pyrolysis_kgps = 0.10,
  .heatRelease_J_per_molFuel = 150e3
});
```

### Safety Invariants (Enforced Throughout Codebase)
- All physical quantities are finite (`std::isfinite()` checks everywhere)
- Non-negative species moles (n_i >= 0)
- Temperature clamped [1.0, 5000.0] K
- All observe() outputs bounded
- No exceptions in production code (all functions `noexcept` or early-return)
- Deterministic RNG (xorshift32 with fixed seed, no `std::random_device`)

### Test Framework Pattern
```cpp
// In TestNumericIntegrity.cpp
static void runNewTest_8X1() {
  // Setup
  Simulation sim;
  sim.resetToScenario({...});
  
  // Execute
  for (int i = 0; i < 100; ++i) {
    sim.step(0.05);
    auto obs = sim.observe();
    
    // Validate
    REQUIRE_FINITE(obs.HRR_W, "8X1 HRR");
    REQUIRE_FINITE(obs.T_K, "8X1 Temperature");
    REQUIRE(obs.T_K >= 1.0 && obs.T_K <= 5000.0, "8X1 T_K in bounds");
  }
  
  std::cout << "[PASS] 8X1 Test description\n";
}

// In main():
runNewTest_8X1();
```

### ValidationSuite Pattern
```cpp
// In ValidationSuite.cpp
void validateNewScenario() {
  Simulation sim;
  sim.resetToScenario({...});
  
  // Run simulation
  sim.commandIgniteOrIncreasePyrolysis();
  for (int i = 0; i < 1200; ++i) {  // 60s at 0.05s timestep
    sim.step(0.05);
  }
  
  auto obs = sim.observe();
  double predicted = obs.T_K;
  double target = 900.0;  // From literature
  double min_lit = 800.0;
  double max_lit = 1000.0;
  
  double error_pct = 100.0 * std::abs(predicted - target) / target;
  bool in_range = (predicted >= min_lit && predicted <= max_lit);
  
  std::cout << "Predicted: " << predicted << " K\n";
  std::cout << "Error: " << error_pct << "%\n";
  std::cout << "Status: " << (in_range ? "PASS" : "FAIL") << "\n";
}
```

---

## Build & Test Commands

```bash
# Build (from d:\Chemsi)
cmake --build build-mingw64 --config Release

# Run tests (from d:\Chemsi\build-mingw64)
.\NumericIntegrity.exe      # Should show 57/57 PASS
.\ValidationSuite.exe        # Should show 4/4 PASS

# After Phase 8 Week 1 complete:
.\NumericIntegrity.exe      # Should show 57/57 PASS (no new tests yet)
.\ValidationSuite.exe        # Should show 7/7 PASS (3 new scenarios)
```

---

## Immediate Next Steps (Week 1, Day 1)

### Task 1: Ship Fire Scenario Implementation
**File to modify:** `cpp_engine/src/ValidationSuite.cpp`

**What to add:**
1. New function `validateShipFire()` following the pattern above
2. Parameters:
   - `volume_m3 = 100.0`
   - `area_m2 = 150.0`
   - `h_W_m2K = 15.0` (aluminum)
   - `ACH = 0.8`
   - `pyrolysis_kgps = 0.10`
   - `heatRelease_J_per_molFuel = 150e3`
3. Run for 60 seconds (1200 steps at dt=0.05s)
4. Check peak temperature against 800-1000K range
5. Calculate error vs midpoint (900K)
6. Output formatted results

**Call site:** Add to `main()` in ValidationSuite.cpp after existing scenarios

**Success criteria:** Ship fire reports <15% error and PASS status

### Task 2: Update ValidationSuite Summary Table
After all 3 scenarios implemented, update the summary table at the end of `main()` to show 7/7 scenarios.

---

## Important Code Patterns to Follow

### 1. Parameter Validation
```cpp
if (!std::isfinite(value) || value < 0.0) {
  return safe_default;  // Never throw
}
```

### 2. Bounds Clamping
```cpp
T_K = std::clamp(T_K, kMinTemp_K, kMaxTemp_K);  // [1.0, 5000.0]
value = std::max(0.0, value);  // Non-negative
```

### 3. Finite Checks in Tests
```cpp
REQUIRE_FINITE(obs.HRR_W, "Test name HRR");
REQUIRE(obs.fuel_kg >= 0.0, "Test name fuel non-negative");
```

### 4. No Exceptions in Production
```cpp
// Good:
if (!isFinitePositive(dt)) return;  // Early return

// Bad:
if (!isFinitePositive(dt)) throw std::runtime_error("Invalid dt");
```

---

## Reference Documents Available

In the `d:\Chemsi` directory:
- `PHASE8_STARTUP.md` - Complete 4-week work plan (545 lines)
- `PHASE8_SESSION_INIT.md` - Session checklist and daily workflow
- `PHASE8_QUICK_REF.md` - One-page Phase 8 summary
- `PHASE7_INTEGRITY_REPORT.md` - Detailed baseline validation report
- `readme.md` - Project overview and current status

---

## Context for AI Assistants

**What I need help with:**
I'm implementing Phase 8 of the VFEP fire simulation project. I need to add 3 new fire scenarios (Ship, Tunnel, Industrial) to expand validation beyond data center fires, then implement a three-zone model and CFD interface.

**Current priority:** Week 1 - Add Ship Fire scenario to ValidationSuite.cpp

**Code style:** 
- C++17, no external dependencies for core physics
- Safety-first (finite checks, bounds clamping, no exceptions)
- Deterministic (fixed RNG seeds, reproducible results)
- Well-commented with physics rationale

**Testing philosophy:**
- Every feature gets numeric integrity tests immediately
- Literature validation required for all scenarios
- No regressions tolerated (all existing tests must continue passing)

**What NOT to do:**
- Do not modify Phase 6 core files (Simulation.cpp, Reactor.cpp, chemistry.cpp, etc.) unless absolutely necessary
- Do not add external dependencies without discussion
- Do not create exceptions in production code
- Do not use `std::random_device` (breaks determinism)

---

## Success Metrics

**Week 1 Complete:**
- ✅ 3 new scenarios implemented in ValidationSuite.cpp
- ✅ All 3 passing with errors <20%
- ✅ ValidationSuite reports 7/7 scenarios PASS
- ✅ 57/57 numeric tests still passing (no regressions)
- ✅ Zero compiler warnings/errors

**Phase 8 Complete (Week 4):**
- ✅ 62/62 numeric tests (5 new for zones + CFD)
- ✅ 7/7 validation scenarios
- ✅ ThreeZoneFireModel operational
- ✅ CFDCoupler interface defined
- ✅ 4 documentation files published

---

## Questions to Ask Me

If you need clarification:
1. "Should I add the scenario to ValidationSuite.cpp or create a new file?"
   → **Add to ValidationSuite.cpp**

2. "What if the scenario doesn't pass immediately?"
   → **Expected. Adjust parameters iteratively within literature bounds.**

3. "Should I add numeric integrity tests now?"
   → **Not yet. Scenarios first (Week 1-2), then zone tests (Week 2-3).**

4. "Can I modify Simulation.cpp?"
   → **Only if adding a new resetToXScenario() helper. Core loop is locked.**

5. "What tolerance is acceptable?"
   → **Ship: <15%, Tunnel: <20%, Industrial: <18% per PHASE8_STARTUP.md**

---

## Ready to Start!

I have:
- ✅ Clean build (57/57 tests, 4/4 scenarios)
- ✅ Phase 8 work plan documented
- ✅ Development environment set up
- ✅ Literature references identified

**First task:** Implement Ship Fire scenario in ValidationSuite.cpp with the parameters specified above.

Please help me add the Ship Fire scenario validation function and integrate it into the ValidationSuite. Let's aim for <15% error vs the 800-1000K literature range (midpoint 900K).
