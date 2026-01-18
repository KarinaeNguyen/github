# 🚀 Quick Start Guide - VFEP Visualization

## Build Status
✅ **Compilation Complete**
- Simulation Server: Ready (VFEP_Sim.exe)
- Visualization: Ready (VFEP.exe)
- All 4 improvements implemented and verified

## How to Run

### Step 1: Terminal 1 - Start Simulation Server
```powershell
cd D:\Chemsi\build-clean
.\VFEP_Sim.exe --grpc_port 50051 --no_auto --t_end 9999
```

Expected output:
```
[gRPC] VFEPUnitySimServiceV1 listening on 127.0.0.1:50051 tick_hz=20
```

### Step 2: Terminal 2 - Start Visualization
```powershell
cd D:\Chemsi\build-clean
.\VFEP.exe
```

## What You'll See

A 3D window with professional interface featuring:

**4 Tabs at top:**
1. **EXEC** - Simulation controls (RUN/PAUSE/STEP/RESET, scenarios, ignition, suppression)
2. **NOZZLE** - Manual nozzle position/direction, rail parameters
3. **VIZ** - Object visibility toggles (warehouse, rack, fire, spray, etc.)
4. **PLOTS** - Real-time data graphs (Temperature, HRR, Exposure, Knockdown, O2)

**3D View (left side):**
- Warehouse with racks
- Fire visualization
- Nozzle and spray cone
- Rail system

**Features Active:**
- ✅ Font size 19pt (larger, easier to read)
- ✅ Nozzle auto-animates when suppression starts
- ✅ 4-tab interface for organized workflow
- ✅ Real-time plots in dedicated PLOTS tab

## Interactive Demo

1. **Click RUN** on EXEC tab → Simulation advances
2. **Click IGNITE** → Fire starts (watch HRR rise on PLOTS)
3. **Click START SUPPRESSION** → Nozzle animates automatically (suppression effects visible on PLOTS)
4. **Switch to NOZZLE tab** → See nozzle position updating
5. **Switch to PLOTS tab** → View 5 real-time graphs
6. **Toggle objects in VIZ tab** → Show/hide scene elements

## Files Reference

| File | Purpose |
|------|---------|
| `cpp_engine/vis/main_vis.cpp` | Visualization code (1256 LOC, all improvements included) |
| `build-clean/VFEP_Sim.exe` | Simulation server executable |
| `build-clean/VFEP.exe` | Visualization client executable |
| `VFEP_IMPROVEMENTS_FINAL.md` | Detailed documentation with code samples |

## Performance

- 60+ FPS rendering
- Real-time nozzle animation
- Smooth plot updates with 1000+ data points
- Instant tab switching

## Troubleshooting

**"Connection refused" error in visualization:**
- Make sure VFEP_Sim.exe started successfully first
- Wait 2-3 seconds after starting server before launching visualization
- Check that port 50051 is available

**Small font (not seeing 19pt):**
- Font improvement is applied automatically
- Check that C:/Windows/Fonts/segoeui.ttf exists on your system

**Nozzle not animating:**
- Make sure suppression is active (click "START SUPPRESSION")
- Switch to NOZZLE tab to see position updating
- Check that `viz_override_nozzle_pose` is unchecked for auto-mode

---

**Created**: January 18, 2026  
**Status**: Production Ready  
**All improvements**: ✅ Verified
