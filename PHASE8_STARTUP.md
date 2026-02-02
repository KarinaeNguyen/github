# Phase 8: Advanced Scenarios & Multi-Zone Physics

**Prepared**: February 2, 2026  
**Status**: READY TO START  
**Foundation**: Phase 7 Week 1 Complete (SensitivityAnalysis, UQ modules, 57/57 tests)

---

## Phase 8 Overview

Phase 8 completes the advanced validation roadmap from Phase 7 by implementing:
1. **New Fire Scenarios** (Ship, Tunnel, Industrial compartments)
2. **True Three-Zone Model** (coupled hot/middle/ambient layers)
3. **CFD Interface** (mock coupling for velocity field import)
4. **Performance Optimization** (<2s per 60s simulation)
5. **Production Documentation** (API reference, user guide, technical report)

### Success Criteria
- ✅ 7+ scenarios validated (currently: 4)
- ✅ Three-zone model operational with validation
- ✅ Performance: <2s per 60s simulation (currently ~0.3s, exceeds target)
- ✅ CFD interface defined with sample import
- ✅ Complete technical documentation suite

---

## What You're Building On

### From Phase 7 (✅ Week 1 Complete)
- **57/57 numeric integrity tests** passing (including 6 Phase 7 tests)
- **4/4 validation scenarios** passing (ISO 9705, NIST, Suppression, Stratification)
- **SensitivityAnalysis module**: Parameter sweeps operational
- **UncertaintyQuantification module**: Monte Carlo + LHS working
- **SweepTool**: CLI utility for batch parameter studies
- Clean build, zero warnings/errors

### Current State
- **NIST Baseline**: 71.4 kW (4.85% error) - LOCKED ✅
- **ISO 9705**: 981K (4.11% error) - LOCKED ✅
- **Suppression**: 79.77% reduction (13.95% error) - LOCKED ✅
- **Stratification**: 272K ΔT (9.26% error) - LOCKED ✅

### Phase 7 Pending Items (Now Phase 8)
- ⏳ New scenarios (Ship Fire, Tunnel Fire, Industrial Compartment)
- ⏳ Three-zone model implementation
- ⏳ CFD coupling interface (mock)
- ⏳ Performance optimization
- ⏳ Extended documentation

---

## Phase 8 Work Plan

### Week 1-2: New Fire Scenarios

#### Scenario 1: Ship Fire (Confined Compartment)
**Objective**: Validate in highly confined space with high heat loss

**Physics**:
- Confined ship compartment (≈100 m³)
- Limited ventilation (0.5-1.0 ACH)
- Aluminum structure (h_W = 15 W/m²K, high conductivity)
- Cargo hold fire (cellulose/paper, 150 kJ/mol)

**Parameters**:
```cpp
// ValidationSuite.cpp: New scenario
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
  
  // Target: Peak T 800-1000K, controlled HRR growth
}
```

**Literature Validation**:
- IMO SOLAS fire testing standards
- Expected peak temperature: 800-1000K
- Time to 500K: <60s
- No flashover (controlled burn)

**Deliverable**: Ship fire scenario passing within ±15% bounds

---

#### Scenario 2: Tunnel Fire (Flow-Driven)
**Objective**: Validate with strong ventilation and stratification

**Physics**:
- Long, slender geometry (500m × 10m² cross-section = 5000 m³)
- Flow-driven ventilation (ACH 2-10, equivalent to 0.5-2.0 m/s)
- Concrete walls (h_W = 8 W/m²K)
- Vehicle fire (hydrocarbons, 450 kJ/mol)

**Parameters**:
```cpp
void validateTunnelFire() {
  Simulation sim;
  sim.resetToScenario({
    .volume_m3 = 5000.0,
    .area_m2 = 1000.0,
    .h_W_m2K = 8.0,  // Concrete
    .ACH = 5.0,      // 1 m/s equivalent
    .pyrolysis_kgps = 0.20,
    .heatRelease_J_per_molFuel = 450e3
  });
  
  // Target: HRR plateau 500-2000 kW, stratified layer >2m
}
```

**Literature Validation**:
- EUREKA 499 tunnel fire experiments
- Memorial Tunnel Fire Ventilation Test Program
- Expected HRR: 500-2000 kW plateau
- Smoke layer height: >2m (good stratification)

**Deliverable**: Tunnel fire scenario passing within ±20% bounds (higher tolerance for flow-dominated case)

---

#### Scenario 3: Industrial Compartment (Warehouse)
**Objective**: Validate large-volume, multi-material fire

**Physics**:
- Large warehouse (2500 m³)
- Mixed ventilation (natural + forced, ACH 1-5)
- Steel frame (h_W = 5 W/m²K, radiates heat)
- Rack storage (mixed materials, 200 kJ/mol)

**Parameters**:
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
  
  // Target: Ceiling T >500K (sprinkler activation)
}
```

**Literature Validation**:
- ISO 9414 (Gross calorific potential)
- ASTM E603 (Room fire test)
- Expected ceiling temperature: >500K (sprinkler trigger)
- Time to untenable conditions: 2-5 min

**Deliverable**: Industrial fire scenario passing within ±18% bounds

---

### Week 2-3: Three-Zone Model

#### Objective
Implement true three-zone fire model with coupled zones (upper hot layer, middle transition, lower ambient).

**Current Limitation**: Single-zone with post-hoc two-zone correction  
**Target**: Coupled three-zone system with mass/energy exchange

#### Architecture

```cpp
// New: ThreeZoneModel.h
namespace vfep::phase8 {

struct Zone {
  double volume_m3;
  double height_m;          // Zone vertical extent
  double T_K;
  double P_Pa;
  std::vector<double> n_mol;  // Species composition per zone
  
  double heatContent_J() const;
  double density_kg_m3() const;
};

class ThreeZoneFireModel {
public:
  ThreeZoneFireModel(double total_height_m, 
                     double floor_area_m2,
                     const std::vector<Species>& species);
  
  void step(double dt, 
            double combustion_HRR_W,
            double cooling_W,
            double ACH);
  
  // Zone access
  const Zone& upperZone() const { return upper_; }
  const Zone& middleZone() const { return middle_; }
  const Zone& lowerZone() const { return lower_; }
  
  // Observables
  double smokeLayerHeight_m() const;
  double averageTemperature_K() const;
  
private:
  Zone upper_;   // Hot smoke layer (0.2-0.4 of height)
  Zone middle_;  // Transition layer (0.3-0.4 of height)
  Zone lower_;   // Cool ambient layer (0.2-0.4 of height)
  
  double total_height_m_;
  double floor_area_m2_;
  
  void updateZoneBoundaries();      // Adjust heights based on density
  void updateMassExchange();        // Buoyancy-driven inter-zone flow
  void updateHeatTransfer();        // Zone-to-zone conduction/radiation
  void updateSpeciesTransport();    // Diffusion across zones
};

} // namespace vfep::phase8
```

#### Implementation Steps

1. **Zone Initialization**
   - Define initial heights: upper (30%), middle (40%), lower (30%)
   - Distribute initial air mass across zones
   - Set initial temperatures (upper=T_amb+10, middle=T_amb+5, lower=T_amb)

2. **Mass Exchange Model**
   ```cpp
   // Buoyancy-driven flow between zones
   double mdot_upper_to_middle = k_exchange * A_interface * sqrt(delta_rho);
   double mdot_middle_to_lower = k_exchange * A_interface * sqrt(delta_rho);
   ```

3. **Heat Transfer**
   ```cpp
   // Conductive/radiative heat exchange
   double Q_upper_to_middle = h_interface * A * (T_upper - T_middle);
   double Q_to_walls = sum_zones(h_W * A_zone * (T_zone - T_amb));
   ```

4. **Species Transport**
   - Concentration-driven diffusion
   - Advective transport with mass flow
   - Smoke filling algorithm

5. **Validation Against ISO 9705**
   - Compare three-zone prediction to single-zone + correction
   - Expected: 1000-1050K (vs 981K current)
   - Better layer height tracking over time

#### Deliverable
- ThreeZoneModel operational
- ISO 9705 validated with three-zone physics
- 3 numeric integrity tests (8A1-8A3) for zone coupling

---

### Week 3-4: CFD Coupling Interface

#### Objective
Define interface for importing CFD velocity fields and exporting VFEP results for comparison.

**Note**: This is a *mock* interface (no actual CFD solver integration). Purpose is to establish data format and comparison methodology.

#### Architecture

```cpp
// New: CFDInterface.h
namespace vfep::phase8 {

struct GridPoint {
  double x, y, z;        // Position (m)
  double T_K;            // Temperature (K)
  double u, v, w;        // Velocity (m/s)
  double rho_kg_m3;      // Density (kg/m³)
};

class CFDCoupler {
public:
  // Import velocity field from CFD output
  void importVelocityField(const std::string& vtk_file);
  void importTemperatureField(const std::string& vtk_file);
  
  // Export VFEP results for comparison
  void exportResults(const std::string& output_vtk);
  void exportComparisonCSV(const std::string& csv_file);
  
  // Query interpolated values at arbitrary points
  double interpolateTemperature(double x, double y, double z, double t);
  Vec3 interpolateVelocity(double x, double y, double z, double t);
  
  // Comparison metrics
  struct ComparisonStats {
    double mean_error;
    double max_error;
    double rmse;
    double correlation;
  };
  
  ComparisonStats compareTemperature(const std::string& cfd_file);
  
private:
  std::vector<GridPoint> grid_;
  
  // Spatial interpolation (trilinear)
  double trilinearInterp(const std::vector<GridPoint>& data,
                         double x, double y, double z);
};

} // namespace vfep::phase8
```

#### File Format Support

**VTK (Visualization Toolkit)** - primary format
```xml
<!-- sample.vtk -->
<VTKFile type="StructuredGrid">
  <StructuredGrid WholeExtent="0 10 0 10 0 10">
    <Piece Extent="0 10 0 10 0 10">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3">
          <!-- x, y, z coordinates -->
        </DataArray>
      </Points>
      <PointData Scalars="Temperature" Vectors="Velocity">
        <DataArray Name="Temperature" type="Float32"/>
        <DataArray Name="Velocity" type="Float32" NumberOfComponents="3"/>
      </PointData>
    </Piece>
  </StructuredGrid>
</VTKFile>
```

**CSV (Simple Format)** - for comparison
```csv
x,y,z,t,T_VFEP,T_CFD,u_CFD,v_CFD,w_CFD
0.5,0.5,0.5,10.0,350.2,348.5,0.1,0.0,0.5
```

#### Mock Comparison Workflow

1. **Generate Mock CFD Data** (for testing)
   ```cpp
   // tools/GenerateMockCFD.cpp
   // Creates synthetic VTK file with analytical temperature/velocity
   ```

2. **Import into VFEP**
   ```cpp
   CFDCoupler coupler;
   coupler.importVelocityField("mock_cfd.vtk");
   ```

3. **Run VFEP Simulation** (with frozen velocity field)

4. **Export and Compare**
   ```cpp
   coupler.exportResults("vfep_results.vtk");
   auto stats = coupler.compareTemperature("mock_cfd.vtk");
   ```

5. **Validation Report**
   - Mean error: ±X K
   - RMSE: Y K
   - Correlation: R² = Z

#### Deliverable
- CFDCoupler class operational
- Mock CFD file generator
- Sample comparison workflow documented
- 2 numeric integrity tests (8B1-8B2) for CFD import/export

---

### Week 4: Performance & Documentation

#### Performance Optimization

**Current Performance**: ~0.3s per 60s simulation (200× real-time) - **Already exceeds <2s target**

**Optional Optimizations** (if needed for batch sweeps):

1. **Vectorize Reaction Kinetics**
   ```cpp
   // Use SIMD for species rate calculations
   #include <immintrin.h>  // AVX/AVX2
   
   // Process 4 species at once with AVX
   __m256d rate = _mm256_mul_pd(k_vec, conc_vec);
   ```

2. **Cache Frequent Calculations**
   ```cpp
   // Cache Arrhenius rates if T_K hasn't changed significantly
   double cached_kT_;
   double cached_T_;
   
   if (abs(T_K - cached_T_) < 1.0) {
     kT = cached_kT_;  // Reuse
   }
   ```

3. **Adaptive Timestepping**
   ```cpp
   // Increase dt during quasi-steady phases
   double adaptiveDt(double base_dt, const Observation& obs) {
     if (obs.HRR_W < 100.0 && obs.T_K < 350.0) {
       return base_dt * 2.0;  // Slow chemistry, larger step OK
     }
     return base_dt;
   }
   ```

**Target**: If needed, optimize batch sweeps (100+ scenarios) to <60s total.

---

#### Documentation Suite

**Phase8_API_Reference.md** - Complete API documentation
```markdown
# VFEP Phase 8 API Reference

## Core Classes
- Simulation
- Reactor
- Chemistry
- Suppression
- Ventilation

## Phase 7 Modules
- SensitivityAnalyzer
- MonteCarloUQ
- SweepTool

## Phase 8 Modules
- ThreeZoneFireModel
- CFDCoupler

## Examples
(Code snippets for common use cases)
```

**Phase8_User_Guide.md** - User workflows
```markdown
# VFEP User Guide

## Quick Start
1. Build and run
2. Run validation scenarios
3. Perform parameter sweeps

## Advanced Workflows
1. Sensitivity analysis
2. Uncertainty quantification
3. Three-zone model usage
4. CFD comparison

## Troubleshooting
```

**Phase8_Technical_Report.md** - Validation report
```markdown
# VFEP Technical Report

## Physics Models
- Combustion chemistry
- Heat transfer
- Ventilation
- Suppression

## Validation Results
- 7 scenarios validated
- Literature comparison
- Uncertainty bounds

## Model Limitations
- Assumptions
- Applicability range
- Future work
```

**Phase8_Scenarios_Catalog.md** - Scenario reference
```markdown
# VFEP Scenarios Catalog

## Scenario 1: ISO 9705 Room Corner Test
- Geometry: 2.4m × 3.6m × 2.4m
- Heat Release: 100 kJ/mol
- Validation: ±5% of 1023K

(Repeat for all 7 scenarios)
```

#### Deliverable
- 4 comprehensive markdown documents
- Code examples in docs/examples/
- Updated readme.md with Phase 8 status

---

## Testing Strategy

### Numeric Integrity Tests (New Phase 8 Tests)

**8A1: Three-Zone Basic Initialization**
```cpp
void runThreeZoneBasicInit() {
  ThreeZoneFireModel tzm(3.0, 10.0, species);
  
  // Verify zone heights sum to total
  REQUIRE(tzm.upperZone().height_m + 
          tzm.middleZone().height_m + 
          tzm.lowerZone().height_m == 3.0);
  
  // Verify initial temperatures finite
  REQUIRE_FINITE(tzm.upperZone().T_K);
  REQUIRE_FINITE(tzm.middleZone().T_K);
  REQUIRE_FINITE(tzm.lowerZone().T_K);
}
```

**8A2: Three-Zone Mass Conservation**
```cpp
void runThreeZoneMassConservation() {
  ThreeZoneFireModel tzm(3.0, 10.0, species);
  
  double total_mass_initial = tzm.totalMass_kg();
  
  for (int i = 0; i < 1000; ++i) {
    tzm.step(0.05, 50e3, 0.0, 2.0);
  }
  
  double total_mass_final = tzm.totalMass_kg();
  
  // With ventilation, mass changes; without, should conserve
  // Test with ACH=0 to verify conservation
}
```

**8A3: Three-Zone Energy Balance**
```cpp
void runThreeZoneEnergyBalance() {
  // Energy in = combustion + inlet
  // Energy out = heat loss + outlet
  // Change in stored energy = Ein - Eout
}
```

**8B1: CFD Import Basic**
```cpp
void runCFDImportBasic() {
  CFDCoupler coupler;
  coupler.importVelocityField("test_data/mock_cfd.vtk");
  
  // Verify grid loaded
  REQUIRE(coupler.gridPointCount() > 0);
  
  // Query temperature at known point
  double T = coupler.interpolateTemperature(0.5, 0.5, 0.5, 10.0);
  REQUIRE_FINITE(T);
}
```

**8B2: CFD Export and Compare**
```cpp
void runCFDExportCompare() {
  CFDCoupler coupler;
  coupler.exportResults("test_output.vtk");
  
  auto stats = coupler.compareTemperature("test_data/reference.vtk");
  REQUIRE_FINITE(stats.mean_error);
  REQUIRE_FINITE(stats.rmse);
}
```

**Target**: 62/62 numeric integrity tests passing (57 current + 5 Phase 8)

---

### Validation Scenarios (Extended Suite)

**New ValidationSuite targets:**
- ✅ ISO 9705 (existing, 4.11% error)
- ✅ NIST Data Center (existing, 4.85% error)
- ✅ Suppression (existing, 13.95% error)
- ✅ Stratification (existing, 9.26% error)
- 🆕 Ship Fire (<15% error vs IMO SOLAS data)
- 🆕 Tunnel Fire (<20% error vs EUREKA 499)
- 🆕 Industrial Fire (<18% error vs ISO 9414)

**Target**: 7/7 scenarios passing

---

## Success Metrics

### Week 1-2 Complete
- ✅ Ship Fire scenario validated (<15% error)
- ✅ Tunnel Fire scenario validated (<20% error)
- ✅ Industrial Fire scenario validated (<18% error)
- ✅ ValidationSuite reports 7/7 scenarios passing

### Week 2-3 Complete
- ✅ ThreeZoneFireModel implemented
- ✅ ISO 9705 validated with three-zone physics (1000-1050K)
- ✅ 3 numeric integrity tests for zones (8A1-8A3)

### Week 3-4 Complete
- ✅ CFDCoupler interface operational
- ✅ Mock CFD import/export working
- ✅ 2 numeric integrity tests for CFD (8B1-8B2)

### Week 4 Complete
- ✅ 4 documentation files published
- ✅ Performance verified (<2s per 60s sim, or optimizations applied)
- ✅ 62/62 numeric integrity tests passing
- ✅ 7/7 validation scenarios passing

---

## Phase 8 Completion Criteria

**All of the following must pass:**

1. ✅ **62/62 numeric integrity tests** (57 existing + 5 Phase 8)
2. ✅ **7/7 validation scenarios** within literature bounds
3. ✅ **Zero compiler warnings/errors**
4. ✅ **ThreeZoneFireModel** operational and validated
5. ✅ **CFDCoupler** interface complete with sample workflow
6. ✅ **Documentation** suite published (4 files)
7. ✅ **Performance** target met (<2s per 60s simulation)

**Phase 8 marks the completion of the advanced validation roadmap and prepares VFEP for production deployment.**

---

## Next Steps (Day 1)

1. **Review Phase 7 status** - confirm all Week 1 items complete
2. **Read this document** end-to-end
3. **Set up development environment** (already done from Phase 7)
4. **Begin Ship Fire scenario implementation** (ValidationSuite.cpp)
5. **Run existing tests** to verify baseline (should be 57/57)

**Let's build Phase 8! 🚀**
