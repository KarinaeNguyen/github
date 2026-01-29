# STL Import Guide for VFEP Visualizer

## Overview
The VFEP visualizer now supports importing 3D models in **Binary STL format**. This allows you to visualize custom room layouts, equipment, racks, or other 3D geometry alongside the fire suppression simulation.

## Supported Format
- **Binary STL** (.stl files)
- ASCII STL is NOT currently supported (binary only)
- File size: Recommended < 10,000 triangles for smooth performance

## How to Use

### 1. Launch the Visualizer
```bash
./VFEP.exe
# or
./VFEP_Vis.exe
```

### 2. Open the STL Tab
In the main UI, click on the **"STL"** tab (between VIZ and PLOTS tabs).

### 3. Load an STL File

**Option A: Manual Path Entry**
1. Type the file path in the "File path" text box
2. Click **"[ LOAD STL ]"** button

**Option B: Quick Load Buttons**
- Click "Load: rack.stl" for a server rack model
- Click "Load: equipment.stl" for equipment cabinet
- Click "Load: room.stl" for room floor outline

### 4. Adjust Transform
- **Position**: Drag the X/Y/Z sliders to move the object
- **Scale**: Adjust the scale slider (0.1 to 10.0)
- Auto-scaling normalizes models to fit the scene

### 5. Display Options
- **Visible**: Toggle object visibility
- **Wireframe**: Switch between solid and wireframe rendering

## Creating STL Files

### Free CAD Software

#### 1. **FreeCAD** (Best for parametric CAD)
- Download: https://www.freecad.org/
- Create your model → File → Export → Select "STL Mesh"
- Choose "Binary" format
- Export to your VFEP directory

#### 2. **Blender** (Best for organic modeling)
- Download: https://www.blender.org/
- Create/import your model
- File → Export → STL (.stl)
- In export settings, ensure "Binary" is selected

#### 3. **TinkerCAD** (Best for beginners, browser-based)
- Website: https://www.tinkercad.com/
- Build your model online
- Export → .STL → Download
- TinkerCAD exports binary STL by default

#### 4. **OpenSCAD** (Best for programmers)
- Download: https://openscad.org/
- Write parametric models in code
- Design → Export as STL

### Python Script Generator
Use the included `create_test_stl.py` script to generate test models:

```bash
python create_test_stl.py
```

This creates:
- `test_cube.stl` - Simple 1m cube
- `rack.stl` - Server rack (0.6 x 2.0 x 0.8m with shelves)
- `equipment.stl` - Equipment cabinet (0.5 x 0.5 x 0.3m)
- `room.stl` - Room floor outline (12 x 6 x 12m)

## Example Use Cases

### Data Center Rack Layout
1. Model your actual server racks in FreeCAD
2. Export as rack.stl
3. Load in VFEP at position matching fire simulation coordinates
4. Visualize fire suppression relative to real equipment

### Custom Room Geometry
1. Create room outline with doors/windows in Blender
2. Export as room.stl
3. Position at origin (0, 0, 0)
4. Scale to match warehouse dimensions

### Equipment Cabinets
1. Model UPS, PDU, or other equipment in TinkerCAD
2. Export and position near fire zones
3. Validate suppression coverage visually

## Technical Details

### Coordinate System
- **X-axis**: Left (-) to Right (+)
- **Y-axis**: Down (-) to Up (+)
- **Z-axis**: Back (-) to Front (+)
- Origin (0, 0, 0) is at warehouse center

### Auto-scaling
- Models are automatically normalized so the largest dimension = 1.0
- Use the Scale slider to adjust relative to scene
- Center point is automatically calculated

### Rendering
- **Solid mode**: Cyan metal material with lighting
- **Wireframe mode**: Bright cyan edges, no lighting
- OpenGL fixed-pipeline rendering for compatibility

### Performance Tips
- Keep triangle count < 10,000 for smooth 60 FPS
- Use wireframe for very complex models
- Binary STL loads ~10x faster than ASCII

## File Paths

### Absolute Paths
```
D:/Chemsi/rack.stl
C:/Users/YourName/Desktop/model.stl
```

### Relative Paths
```
rack.stl                    # Same directory as VFEP.exe
../models/equipment.stl     # Parent directory
./stl/room.stl              # Subdirectory
```

## Troubleshooting

### "Failed to open STL file"
- Check file path is correct (use absolute path)
- Ensure file exists and has .stl extension
- Verify you have read permissions

### "Invalid triangle count"
- File may be ASCII STL (not supported yet)
- File may be corrupted
- Try re-exporting from your CAD software

### Model appears too large/small
- Adjust the **Scale** slider in the STL tab
- Check your CAD software units (meters recommended)

### Model not visible
- Click "Visible" checkbox in STL tab
- Check if model is inside other geometry
- Try adjusting Position sliders

### Model looks wrong
- Try wireframe mode to check geometry
- Verify STL export settings in your CAD software
- Ensure binary format (not ASCII)

## Known Limitations
- Binary STL only (ASCII not yet supported)
- One mesh at a time (no multi-object support yet)
- No texture/color import (uses default material)
- No animation support
- Fixed-pipeline rendering (no advanced shaders)

## Planned Features
- Multiple STL objects simultaneously
- ASCII STL support
- OBJ file format support
- Material/color preservation
- Drag-and-drop file loading
- STL file browser dialog

## Example Workflow

### Creating a Custom Rack Model

1. **Measure your real rack**: 600mm (W) × 2000mm (H) × 800mm (D)

2. **Open FreeCAD**:
   - Create → Part → Box
   - Set dimensions: 0.6 × 2.0 × 0.8 meters
   - Add shelves using Part → Box for each shelf
   - Union all parts

3. **Export**:
   - File → Export
   - Type: "STL Mesh (*.stl)"
   - Format: Binary
   - Save as: server_rack.stl

4. **Load in VFEP**:
   - Open VFEP visualizer
   - Go to STL tab
   - Enter path: D:/Chemsi/server_rack.stl
   - Click [ LOAD STL ]

5. **Position**:
   - Set Position: (0, 1, 0) to center at rack height
   - Adjust scale if needed
   - Toggle wireframe to check alignment

6. **Verify Coverage**:
   - Run simulation
   - Watch spray cone interaction with imported rack
   - Validate suppression effectiveness

## SolidWorks Files

**Q: Can I use SolidWorks (.sldprt) files?**

**A: Not directly.** SolidWorks files must be converted to STL first:

### SolidWorks Export Instructions:
1. Open your part/assembly in SolidWorks
2. File → Save As
3. Save as type: **"STL (*.stl)"**
4. Options → Output as: **Binary**
5. Units: **Meters** (recommended)
6. Resolution: **Fine** or **Custom**
7. Save

Then load the .stl file in VFEP as described above.

### Alternative: Use Free Converters
If you don't have SolidWorks:
- **FreeCAD** can import some SolidWorks files
- **Online converters**: Several websites convert SLDPRT → STL
- Ask someone with SolidWorks to export for you

## Support

For issues or feature requests, check:
- Build logs in the cpp_engine/vis directory
- Console output when loading STL files
- Verify OpenGL version supports required features

---
**Last Updated**: January 29, 2026  
**Version**: 1.0 (Initial STL support)
