# Phase 8 Quick Reference

**Phase**: Advanced Scenarios & Multi-Zone Physics  
**Duration**: 4 weeks  
**Foundation**: Phase 7 Week 1 complete (57/57 tests, 4/4 scenarios)

---

## Quick Goals

✅ **New Scenarios**: Ship, Tunnel, Industrial (3 scenarios → 7 total)  
✅ **Three-Zone Model**: Coupled hot/middle/ambient layers  
✅ **CFD Interface**: Mock velocity field import/export  
✅ **Documentation**: 4 complete reference documents  

---

## Week-by-Week

### Week 1-2: New Scenarios
- Ship Fire (confined, aluminum, 800-1000K)
- Tunnel Fire (flow-driven, 500-2000 kW)
- Industrial Fire (large volume, >500K ceiling)

### Week 2-3: Three-Zone Model
- Zone coupling (mass/energy exchange)
- ISO 9705 validation (1000-1050K)
- 3 numeric tests (8A1-8A3)

### Week 3-4: CFD & Docs
- VTK import/export
- Mock CFD workflow
- 4 documentation files
- 2 numeric tests (8B1-8B2)

---

## Target Metrics

**Tests**: 62/62 (57 current + 5 Phase 8)  
**Scenarios**: 7/7 (4 current + 3 Phase 8)  
**Performance**: <2s per 60s sim (currently ~0.3s ✅)  
**Docs**: 4 files complete  

---

## New Files to Create

```
cpp_engine/include/ThreeZoneModel.h
cpp_engine/src/ThreeZoneModel.cpp
cpp_engine/include/CFDInterface.h
cpp_engine/src/CFDInterface.cpp
cpp_engine/tools/GenerateMockCFD.cpp

PHASE8_API_Reference.md
PHASE8_User_Guide.md
PHASE8_Technical_Report.md
PHASE8_Scenarios_Catalog.md
```

---

## First Steps

1. **Verify baseline**: Run `.\NumericIntegrity.exe` (57/57)
2. **Start Ship Fire**: Add to ValidationSuite.cpp
3. **Parameters**:
   - Volume: 100 m³
   - Area: 150 m²
   - h_W: 15 W/m²K (aluminum)
   - ACH: 0.8
   - Pyrolysis: 0.10 kg/s
   - Heat Release: 150 kJ/mol
4. **Target**: 800-1000K peak temperature (<15% error)

---

## Quick Commands

```bash
# Build
cmake --build build-mingw64 --config Release

# Test
cd build-mingw64
.\NumericIntegrity.exe    # 57/57 → 62/62
.\ValidationSuite.exe      # 4/4 → 7/7
```

---

## Literature Sources

- **Ship Fire**: IMO SOLAS fire testing
- **Tunnel Fire**: EUREKA 499, Memorial Tunnel
- **Industrial**: ISO 9414, ASTM E603

---

## Success = Phase 8 Complete

✅ 62/62 tests  
✅ 7/7 scenarios  
✅ Three-zone operational  
✅ CFD interface defined  
✅ 4 docs published  
✅ Zero warnings/errors  

**Ready to start!** 🚀
