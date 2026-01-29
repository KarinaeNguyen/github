# VFEP STL Import - Quick Reference

## ✅ YES - This Works!

### Supported Formats
- ✅ **Binary STL** (.stl) - FULLY SUPPORTED
- ❌ ASCII STL - Not yet supported
- ❌ SolidWorks (.sldprt) - Must convert to STL first

### Best Free Alternatives to SolidWorks

| Software | Type | Best For | Export to STL |
|----------|------|----------|---------------|
| **FreeCAD** | Parametric CAD | Engineering/Mechanical | ✅ Yes (Binary) |
| **Blender** | 3D Modeling | Organic/Complex shapes | ✅ Yes (Binary) |
| **TinkerCAD** | Browser CAD | Beginners/Quick models | ✅ Yes (Binary) |
| **OpenSCAD** | Code-based CAD | Programmers | ✅ Yes (Binary) |

## 🚀 Quick Start (3 Steps)

1. **Get an STL file**
   ```bash
   python create_test_stl.py  # Creates example files
   ```

2. **Launch VFEP**
   ```bash
   ./VFEP.exe
   ```

3. **Load in UI**
   - Click **"STL"** tab
   - Enter file path: `rack.stl`
   - Click **[ LOAD STL ]**

## 💡 Use Cases

- ✅ Import actual data center rack models
- ✅ Visualize room layouts with equipment
- ✅ Validate fire suppression coverage
- ✅ Test different equipment configurations
- ✅ Create custom scenarios

## ⚙️ Controls (STL Tab)

```
File path: [____________]  <- Type path here
           [ LOAD STL ]    <- Click to load

Transform:
  Position: [X] [Y] [Z]    <- Move object
  Scale:    [0.1 - 10.0]   <- Resize

Display:
  [✓] Visible              <- Show/hide
  [ ] Wireframe            <- Solid or edges only
```

## 📝 Example Paths

Windows:
```
D:\Chemsi\rack.stl
C:\Models\equipment.stl
```

Relative:
```
rack.stl
../stl/room.stl
```

## 🔧 Convert SolidWorks to STL

**In SolidWorks:**
1. File → Save As
2. Type: **STL (*.stl)**
3. Options → **Binary**
4. Units: **Meters**
5. Save

**Without SolidWorks:**
- Use FreeCAD to open .sldprt
- Export as STL

## 📦 Included Test Files

Run `python create_test_stl.py` to generate:

- `test_cube.stl` - 1×1×1 meter cube
- `rack.stl` - Server rack with shelves (0.6×2.0×0.8m)
- `equipment.stl` - Small cabinet (0.5×0.5×0.3m)
- `room.stl` - Floor outline (12×6×12m)

## ⚡ Performance

| Triangle Count | Performance |
|----------------|-------------|
| < 1,000 | Excellent (60+ FPS) |
| 1,000 - 10,000 | Good (30-60 FPS) |
| > 10,000 | Use wireframe mode |

## 🎨 Rendering

- **Solid**: Cyan metal material with lighting
- **Wireframe**: Bright cyan edges
- Auto-centers and auto-scales models

## ❓ Common Issues

**"File not found"**
→ Use absolute path: `D:\Chemsi\rack.stl`

**"Invalid triangle count"**
→ File is ASCII STL (not supported) or corrupted
→ Re-export as Binary STL

**Model too small/large**
→ Adjust Scale slider

**Model not visible**
→ Check "Visible" checkbox
→ Adjust Position sliders

## 📚 Full Documentation

See **STL_IMPORT_GUIDE.md** for complete details.

---

**Summary:** Yes, you can import 3D sketches! Use **Binary STL format** from any free CAD software. SolidWorks files need conversion to STL first.
