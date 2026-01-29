#!/usr/bin/env python3
"""
Simple script to create binary STL files for testing the VFEP visualizer.
Creates basic geometric shapes in STL format.
"""

import struct
import math

def write_binary_stl(filename, triangles):
    """
    Write triangles to binary STL format.
    triangles: list of (normal, v0, v1, v2) where each is (x,y,z)
    """
    with open(filename, 'wb') as f:
        # 80-byte header
        f.write(b'Binary STL created by Python' + b'\x00' * 52)
        
        # Number of triangles
        f.write(struct.pack('<I', len(triangles)))
        
        # Write each triangle
        for normal, v0, v1, v2 in triangles:
            # Normal (3 floats)
            f.write(struct.pack('<fff', *normal))
            # Vertex 0 (3 floats)
            f.write(struct.pack('<fff', *v0))
            # Vertex 1 (3 floats)
            f.write(struct.pack('<fff', *v1))
            # Vertex 2 (3 floats)
            f.write(struct.pack('<fff', *v2))
            # Attribute byte count (uint16)
            f.write(struct.pack('<H', 0))

def compute_normal(v0, v1, v2):
    """Compute triangle normal using cross product"""
    ax = v1[0] - v0[0]
    ay = v1[1] - v0[1]
    az = v1[2] - v0[2]
    
    bx = v2[0] - v0[0]
    by = v2[1] - v0[1]
    bz = v2[2] - v0[2]
    
    nx = ay * bz - az * by
    ny = az * bx - ax * bz
    nz = ax * by - ay * bx
    
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length > 0:
        nx, ny, nz = nx/length, ny/length, nz/length
    
    return (nx, ny, nz)

def create_cube_stl(filename, size=1.0):
    """Create a cube STL file"""
    s = size / 2
    
    # 8 vertices of a cube
    v = [
        (-s, -s, -s), (s, -s, -s), (s, s, -s), (-s, s, -s),  # back face
        (-s, -s, s), (s, -s, s), (s, s, s), (-s, s, s)       # front face
    ]
    
    # 12 triangles (2 per face, 6 faces)
    faces = [
        # Back face (z = -s)
        (v[0], v[1], v[2]), (v[0], v[2], v[3]),
        # Front face (z = s)
        (v[4], v[6], v[5]), (v[4], v[7], v[6]),
        # Left face (x = -s)
        (v[0], v[3], v[7]), (v[0], v[7], v[4]),
        # Right face (x = s)
        (v[1], v[5], v[6]), (v[1], v[6], v[2]),
        # Bottom face (y = -s)
        (v[0], v[4], v[5]), (v[0], v[5], v[1]),
        # Top face (y = s)
        (v[3], v[2], v[6]), (v[3], v[6], v[7]),
    ]
    
    triangles = []
    for v0, v1, v2 in faces:
        normal = compute_normal(v0, v1, v2)
        triangles.append((normal, v0, v1, v2))
    
    write_binary_stl(filename, triangles)
    print(f"Created {filename} with {len(triangles)} triangles")

def create_rack_stl(filename):
    """Create a simple server rack shape"""
    triangles = []
    
    # Rack dimensions (meters): 0.6 wide, 2.0 tall, 0.8 deep
    w, h, d = 0.6, 2.0, 0.8
    
    # Create outer box
    hw, hh, hd = w/2, h/2, d/2
    v = [
        (-hw, -hh, -hd), (hw, -hh, -hd), (hw, hh, -hd), (-hw, hh, -hd),
        (-hw, -hh, hd), (hw, -hh, hd), (hw, hh, hd), (-hw, hh, hd)
    ]
    
    faces = [
        # Back, Front, Left, Right, Bottom, Top
        (v[0], v[1], v[2]), (v[0], v[2], v[3]),
        (v[4], v[6], v[5]), (v[4], v[7], v[6]),
        (v[0], v[3], v[7]), (v[0], v[7], v[4]),
        (v[1], v[5], v[6]), (v[1], v[6], v[2]),
        (v[0], v[4], v[5]), (v[0], v[5], v[1]),
        (v[3], v[2], v[6]), (v[3], v[6], v[7]),
    ]
    
    for v0, v1, v2 in faces:
        normal = compute_normal(v0, v1, v2)
        triangles.append((normal, v0, v1, v2))
    
    # Add horizontal shelves (5 shelves)
    for i in range(5):
        y = -hh + (i + 1) * h / 6
        thickness = 0.02
        
        shelf_verts = [
            (-hw, y - thickness, -hd), (hw, y - thickness, -hd),
            (hw, y + thickness, -hd), (-hw, y + thickness, -hd),
            (-hw, y - thickness, hd), (hw, y - thickness, hd),
            (hw, y + thickness, hd), (-hw, y + thickness, hd)
        ]
        
        shelf_faces = [
            (shelf_verts[0], shelf_verts[1], shelf_verts[2]), (shelf_verts[0], shelf_verts[2], shelf_verts[3]),
            (shelf_verts[4], shelf_verts[6], shelf_verts[5]), (shelf_verts[4], shelf_verts[7], shelf_verts[6]),
            (shelf_verts[0], shelf_verts[3], shelf_verts[7]), (shelf_verts[0], shelf_verts[7], shelf_verts[4]),
            (shelf_verts[1], shelf_verts[5], shelf_verts[6]), (shelf_verts[1], shelf_verts[6], shelf_verts[2]),
        ]
        
        for v0, v1, v2 in shelf_faces:
            normal = compute_normal(v0, v1, v2)
            triangles.append((normal, v0, v1, v2))
    
    write_binary_stl(filename, triangles)
    print(f"Created {filename} with {len(triangles)} triangles")

def create_equipment_stl(filename):
    """Create a simple equipment cabinet"""
    triangles = []
    
    # Equipment dimensions: 0.5 x 0.5 x 0.3
    w, h, d = 0.5, 0.5, 0.3
    hw, hh, hd = w/2, h/2, d/2
    
    v = [
        (-hw, -hh, -hd), (hw, -hh, -hd), (hw, hh, -hd), (-hw, hh, -hd),
        (-hw, -hh, hd), (hw, -hh, hd), (hw, hh, hd), (-hw, hh, hd)
    ]
    
    faces = [
        (v[0], v[1], v[2]), (v[0], v[2], v[3]),
        (v[4], v[6], v[5]), (v[4], v[7], v[6]),
        (v[0], v[3], v[7]), (v[0], v[7], v[4]),
        (v[1], v[5], v[6]), (v[1], v[6], v[2]),
        (v[0], v[4], v[5]), (v[0], v[5], v[1]),
        (v[3], v[2], v[6]), (v[3], v[6], v[7]),
    ]
    
    for v0, v1, v2 in faces:
        normal = compute_normal(v0, v1, v2)
        triangles.append((normal, v0, v1, v2))
    
    write_binary_stl(filename, triangles)
    print(f"Created {filename} with {len(triangles)} triangles")

def create_room_stl(filename):
    """Create a simple room outline"""
    triangles = []
    
    # Room dimensions: 12 x 6 x 12 (matching warehouse_half * 2)
    w, h, d = 12.0, 6.0, 12.0
    hw, hh, hd = w/2, h/2, d/2
    thickness = 0.1
    
    # Floor
    floor_verts = [
        (-hw, -hh, -hd), (hw, -hh, -hd), (hw, -hh, hd), (-hw, -hh, hd),
        (-hw, -hh + thickness, -hd), (hw, -hh + thickness, -hd),
        (hw, -hh + thickness, hd), (-hw, -hh + thickness, hd)
    ]
    
    floor_faces = [
        (floor_verts[0], floor_verts[2], floor_verts[1]), (floor_verts[0], floor_verts[3], floor_verts[2]),
        (floor_verts[4], floor_verts[5], floor_verts[6]), (floor_verts[4], floor_verts[6], floor_verts[7]),
    ]
    
    for v0, v1, v2 in floor_faces:
        normal = compute_normal(v0, v1, v2)
        triangles.append((normal, v0, v1, v2))
    
    write_binary_stl(filename, triangles)
    print(f"Created {filename} with {len(triangles)} triangles")

if __name__ == '__main__':
    create_cube_stl('test_cube.stl', 1.0)
    create_rack_stl('rack.stl')
    create_equipment_stl('equipment.stl')
    create_room_stl('room.stl')
    print("\nAll test STL files created successfully!")
    print("Copy these files to your VFEP build directory to test.")
