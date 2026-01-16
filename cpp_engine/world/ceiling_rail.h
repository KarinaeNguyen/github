#pragma once

// world/ceiling_rail.h
//
// Step 1 (architecture): First-class, room-attached ceiling rail geometry.
//
// Design goals:
//   - No ImGui / ImPlot / OpenGL dependencies.
//   - Deterministic geometry driven by room + rack inputs.
//   - Simple interface for visualization to render the rail (corners/segments).
//   - Future-ready: exposes a perimeter parameterization (s along the loop)
//     and a stable XZ projection primitive for control/kinematics.
//
// Coordinate convention (matches current main_vis.cpp):
//   - Y is up.
//   - Rail is an axis-aligned rectangle in the XZ plane.
//   - Default ceiling height is derived from a "warehouse half" box where
//       floor_y = 0 and ceiling_y = 2 * warehouse_half.y.
//     Alternatively, you may pass an explicit ceiling_y_m in inputs.
//
// Corner ordering (around the perimeter):
//   0: (x0, z0)
//   1: (x1, z0)
//   2: (x1, z1)
//   3: (x0, z1)
// Segments:
//   seg 0: 0 -> 1 (+X)
//   seg 1: 1 -> 2 (+Z)
//   seg 2: 2 -> 3 (-X)
//   seg 3: 3 -> 0 (-Z)

#include <array>
#include <cstddef>

namespace vfep {
namespace world {

struct Vec3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct CeilingRailConfig {
    // Vertical offset from the room ceiling. Default = 0.30 m (30 cm).
    double drop_from_ceiling_m = 0.30;

    // Extra clearance beyond the rack footprint in X and Z.
    // Default matches current main_vis.cpp.
    double margin_from_rack_m = 0.25;
};

struct CeilingRailInputs {
    // If ceiling_y_m > 0, it is used directly.
    // Otherwise, ceiling_y is derived as: ceiling_y = 2 * warehouse_half_m.y.
    double ceiling_y_m = 0.0;

    // Room visualization convenience. Only warehouse_half_m.y is required if
    // ceiling_y_m == 0; other components are not used by the rail.
    Vec3d warehouse_half_m {6.0, 3.0, 6.0};

    // Rack footprint that drives the rail extents.
    Vec3d rack_center_m {0.0, 1.0, 0.0};
    Vec3d rack_half_m   {0.6, 1.0, 0.4};
};

struct CeilingRailGeometry {
    static constexpr int kNumCorners  = 4;
    static constexpr int kNumSegments = 4;

    // Rail height in room coordinates.
    double y_m = 0.0;

    // Axis-aligned extents in the XZ plane.
    double x0_m = 0.0;
    double x1_m = 0.0;
    double z0_m = 0.0;
    double z1_m = 0.0;

    // Corners in room coordinates (see ordering in header comment).
    std::array<Vec3d, kNumCorners> corners_room_m {};

    // Segment lengths for the ordered loop.
    std::array<double, kNumSegments> seg_len_m {{0.0, 0.0, 0.0, 0.0}};

    // Total perimeter length.
    double perimeter_m = 0.0;

    // Convenience: returns corners as 4 segments (a,b) for rendering.
    std::array<std::array<Vec3d, 2>, kNumSegments> segments_room_m() const {
        return { {
            {corners_room_m[0], corners_room_m[1]},
            {corners_room_m[1], corners_room_m[2]},
            {corners_room_m[2], corners_room_m[3]},
            {corners_room_m[3], corners_room_m[0]},
        } };
    }
};

class CeilingRail {
public:
    struct Projection {
        // Arc-length parameter along the perimeter in [0, perimeter).
        double s_m = 0.0;

        // Closest point on the rail (room coordinates).
        Vec3d pos_room_m {};

        // Euclidean distance in XZ from query point to rail (>= 0).
        double dist_m = 0.0;

        // Segment index where the closest point lies (0..3).
        int segment_idx = 0;
    };

    CeilingRail() = default;
    explicit CeilingRail(const CeilingRailConfig& cfg) : cfg_(cfg) {}

    void setConfig(const CeilingRailConfig& cfg) { cfg_ = cfg; }
    const CeilingRailConfig& config() const { return cfg_; }

    // Recompute the rail geometry from inputs.
    // This is intentionally side-effect free beyond updating geometry().
    void recompute(const CeilingRailInputs& in);

    bool isValid() const { return valid_; }
    const CeilingRailGeometry& geometry() const { return geo_; }

    // Wrap s into [0, perimeter). If the rail is invalid, returns 0.
    double wrapS(double s_m) const;

    // Signed smallest delta from s_from to s_to along the loop.
    // Range is approximately (-perimeter/2, +perimeter/2].
    // If the rail is invalid, returns 0.
    double shortestDeltaS(double s_from_m, double s_to_m) const;

    // Evaluate position on rail at arc-length s.
    // If invalid, returns {0,0,0}.
    Vec3d evalPosition(double s_m) const;

    // Evaluate tangent direction (unit-length, y=0) at arc-length s.
    // If invalid, returns {1,0,0}.
    Vec3d evalTangent(double s_m) const;

    // Project an XZ point onto the rail, returning the nearest point.
    // Optional s_hint_m is used for deterministic tie-breaking near corners.
    Projection projectNearestXZ(double x_m, double z_m, double s_hint_m = 0.0) const;

private:
    CeilingRailConfig cfg_ {};
    CeilingRailGeometry geo_ {};
    bool valid_ = false;
};

} // namespace world
} // namespace vfep

