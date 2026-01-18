// world/ceiling_rail.cpp
//
// Step 1 (architecture): First-class, room-attached ceiling rail geometry.
//
// Implementation notes:
//   - Matches the prior UI-only rail behavior in main_vis.cpp:
//       * ceiling_y = (inputs.ceiling_y_m > 0) ? inputs.ceiling_y_m : 2 * warehouse_half_m.y
//       * rail_y    = ceiling_y - cfg.drop_from_ceiling_m
//       * x extents follow rack footprint +/- (rack_half.x + margin)
//       * z extents follow rack footprint +/- (rack_half.z + margin)
//   - No dependencies on ImGui / ImPlot / OpenGL.
//   - Deterministic and side-effect free beyond updating internal geometry.

#include "ceiling_rail.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace vfep {
namespace world {

// Numeric stability thresholds
static constexpr double kEpsilon = 1e-12;
static constexpr double kMinGeometry = 1e-6;  // Minimum dimension before considering degenerate
static constexpr double kMinPerimeter = 1e-12; // Minimum perimeter for valid rail

static inline double clampd(double v, double lo, double hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

static inline double sqr(double v) { return v * v; }

static inline Vec3d make_v3(double x, double y, double z) {
    Vec3d out;
    out.x = x;
    out.y = y;
    out.z = z;
    return out;
}

static inline Vec3d sub_v3(const Vec3d& a, const Vec3d& b) {
    return make_v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline Vec3d norm_xz(const Vec3d& v, const Vec3d& fallback) {
    const double l2 = v.x * v.x + v.z * v.z;
    if (l2 <= kEpsilon * kEpsilon) {
        return fallback;
    }
    const double inv_l = 1.0 / std::sqrt(l2);
    return make_v3(v.x * inv_l, 0.0, v.z * inv_l);
}

void CeilingRail::recompute(const CeilingRailInputs& in) {
    // Input validation
    valid_ = false;
    
    // Validate rack has positive dimensions
    if (in.rack_half_m.x <= kMinGeometry || in.rack_half_m.z <= kMinGeometry) {
        return;  // Degenerate rack
    }
    
    // Validate warehouse has positive height
    if (in.warehouse_half_m.y <= kMinGeometry) {
        return;  // Degenerate warehouse
    }
    
    // Ceiling height must be positive
    double ceiling_y = (in.ceiling_y_m > 0.0)
        ? in.ceiling_y_m
        : (2.0 * in.warehouse_half_m.y);
    
    if (ceiling_y <= 0.0) {
        return;  // Invalid ceiling height
    }
    
    // Rail height is ceiling minus drop (should be positive)
    geo_.y_m = ceiling_y - cfg_.drop_from_ceiling_m;
    
    // Allow rail slightly below floor for visualization, but reject extreme cases
    if (std::abs(geo_.y_m) > 1000.0) {
        return;  // Unreasonable geometry
    }

    // 3) Extents driven by rack footprint (matches main_vis.cpp).
    geo_.x0_m = in.rack_center_m.x - (in.rack_half_m.x + cfg_.margin_from_rack_m);
    geo_.x1_m = in.rack_center_m.x + (in.rack_half_m.x + cfg_.margin_from_rack_m);
    geo_.z0_m = in.rack_center_m.z - (in.rack_half_m.z + cfg_.margin_from_rack_m);
    geo_.z1_m = in.rack_center_m.z + (in.rack_half_m.z + cfg_.margin_from_rack_m);

    // Ensure canonical ordering (x0<=x1, z0<=z1).
    if (geo_.x0_m > geo_.x1_m) std::swap(geo_.x0_m, geo_.x1_m);
    if (geo_.z0_m > geo_.z1_m) std::swap(geo_.z0_m, geo_.z1_m);

    // 4) Corners (ordering documented in the header).
    geo_.corners_room_m[0] = make_v3(geo_.x0_m, geo_.y_m, geo_.z0_m);
    geo_.corners_room_m[1] = make_v3(geo_.x1_m, geo_.y_m, geo_.z0_m);
    geo_.corners_room_m[2] = make_v3(geo_.x1_m, geo_.y_m, geo_.z1_m);
    geo_.corners_room_m[3] = make_v3(geo_.x0_m, geo_.y_m, geo_.z1_m);

    // 5) Segment lengths and perimeter.
    const double len_x = std::abs(geo_.x1_m - geo_.x0_m);
    const double len_z = std::abs(geo_.z1_m - geo_.z0_m);

    geo_.seg_len_m[0] = len_x;
    geo_.seg_len_m[1] = len_z;
    geo_.seg_len_m[2] = len_x;
    geo_.seg_len_m[3] = len_z;

    geo_.perimeter_m = 2.0 * (len_x + len_z);

    // Valid only if we have meaningful perimeter.
    valid_ = (geo_.perimeter_m > kMinPerimeter);
}

double CeilingRail::wrapS(double s_m) const {
    if (!valid_ || geo_.perimeter_m <= 0.0) {
        return 0.0;
    }
    const double p = geo_.perimeter_m;
    double w = std::fmod(s_m, p);
    if (w < 0.0) w += p;
    // Guard against fmod returning p due to numeric edge cases.
    if (w >= p - kEpsilon) w = 0.0;
    return w;
}

double CeilingRail::shortestDeltaS(double s_from_m, double s_to_m) const {
    if (!valid_ || geo_.perimeter_m <= 0.0) {
        return 0.0;
    }
    const double p = geo_.perimeter_m;
    const double a = wrapS(s_from_m);
    const double b = wrapS(s_to_m);
    double d = b - a;

    const double half = 0.5 * p;
    if (d > half) {
        d -= p;
    } else if (d <= -half) {
        d += p;
    }
    return d;
}

// Internal: locate segment index for s (wrapped), and param u in [0,1] along that segment.
static inline int segment_from_s(const CeilingRailGeometry& g, double s_wrapped, double& u01) {
    // We assume g.perimeter_m > 0 and s_wrapped in [0, perimeter).
    double accum = 0.0;
    for (int seg = 0; seg < CeilingRailGeometry::kNumSegments; ++seg) {
        const double L = g.seg_len_m[seg];
        const double next = accum + L;

        // Use <= so exact corner values map deterministically to the earlier segment.
        // For the last segment, always select it if not selected earlier.
        const bool is_last = (seg == CeilingRailGeometry::kNumSegments - 1);
        if (is_last || s_wrapped <= next) {
            if (L > kEpsilon) {
                u01 = (s_wrapped - accum) / L;
                u01 = clampd(u01, 0.0, 1.0);
            } else {
                u01 = 0.0;
            }
            return seg;
        }
        accum = next;
    }

    // Fallback (should never happen).
    u01 = 0.0;
    return 0;
}

Vec3d CeilingRail::evalPosition(double s_m) const {
    if (!valid_) {
        return make_v3(0.0, 0.0, 0.0);
    }

    const double s = wrapS(s_m);
    double u = 0.0;
    const int seg = segment_from_s(geo_, s, u);

    // Segment endpoints by index.
    const int i0 = seg;
    const int i1 = (seg + 1) % CeilingRailGeometry::kNumCorners;

    const Vec3d& a = geo_.corners_room_m[(std::size_t)i0];
    const Vec3d& b = geo_.corners_room_m[(std::size_t)i1];

    const double x = a.x + u * (b.x - a.x);
    const double z = a.z + u * (b.z - a.z);
    return make_v3(x, geo_.y_m, z);
}

Vec3d CeilingRail::evalTangent(double s_m) const {
    if (!valid_) {
        return make_v3(1.0, 0.0, 0.0);
    }

    const double s = wrapS(s_m);
    double u = 0.0;
    const int seg = segment_from_s(geo_, s, u);

    (void)u; // u not needed for tangent in the axis-aligned case.

    const int i0 = seg;
    const int i1 = (seg + 1) % CeilingRailGeometry::kNumCorners;

    const Vec3d& a = geo_.corners_room_m[(std::size_t)i0];
    const Vec3d& b = geo_.corners_room_m[(std::size_t)i1];

    const Vec3d d = sub_v3(b, a);
    return norm_xz(d, make_v3(1.0, 0.0, 0.0));
}

CeilingRail::Projection CeilingRail::projectNearestXZ(double x_m, double z_m, double s_hint_m) const {
    Projection out;
    if (!valid_) {
        out.s_m = 0.0;
        out.pos_room_m = make_v3(0.0, 0.0, 0.0);
        out.dist_m = 0.0;
        out.segment_idx = 0;
        return out;
    }

    struct Cand {
        double s_m_raw = 0.0; // raw in [0, perimeter] (may equal perimeter at closing corner)
        Vec3d pos_room_m {};
        double d2 = 0.0;
        int seg = 0;
    };

    const double x0 = geo_.x0_m;
    const double x1 = geo_.x1_m;
    const double z0 = geo_.z0_m;
    const double z1 = geo_.z1_m;

    const double len0 = geo_.seg_len_m[0];
    const double len1 = geo_.seg_len_m[1];
    const double len2 = geo_.seg_len_m[2];
    const double len3 = geo_.seg_len_m[3];

    const double dx = x1 - x0;
    const double dz = z1 - z0;

    Cand cands[CeilingRailGeometry::kNumSegments];

    // --- Segment 0: 0 -> 1 (z = z0, x in [x0,x1])
    {
        const double xc = clampd(x_m, x0, x1);
        const double zc = z0;
        const double denom = dx;
        double u = 0.0;
        if (std::abs(denom) > 1e-24) {
            u = (xc - x0) / denom;
            u = clampd(u, 0.0, 1.0);
        }
        cands[0].s_m_raw = u * len0;
        cands[0].pos_room_m = make_v3(xc, geo_.y_m, zc);
        cands[0].d2 = sqr(x_m - xc) + sqr(z_m - zc);
        cands[0].seg = 0;
    }

    // --- Segment 1: 1 -> 2 (x = x1, z in [z0,z1])
    {
        const double xc = x1;
        const double zc = clampd(z_m, z0, z1);
        const double denom = dz;
        double u = 0.0;
        if (std::abs(denom) > 1e-24) {
            u = (zc - z0) / denom;
            u = clampd(u, 0.0, 1.0);
        }
        cands[1].s_m_raw = len0 + u * len1;
        cands[1].pos_room_m = make_v3(xc, geo_.y_m, zc);
        cands[1].d2 = sqr(x_m - xc) + sqr(z_m - zc);
        cands[1].seg = 1;
    }

    // --- Segment 2: 2 -> 3 (z = z1, x in [x0,x1], traveling toward -X)
    {
        const double xc = clampd(x_m, x0, x1);
        const double zc = z1;
        const double denom = dx;
        double u = 0.0;
        if (std::abs(denom) > 1e-24) {
            // u=0 at x=x1, u=1 at x=x0
            u = (x1 - xc) / denom;
            u = clampd(u, 0.0, 1.0);
        }
        cands[2].s_m_raw = len0 + len1 + u * len2;
        cands[2].pos_room_m = make_v3(xc, geo_.y_m, zc);
        cands[2].d2 = sqr(x_m - xc) + sqr(z_m - zc);
        cands[2].seg = 2;
    }

    // --- Segment 3: 3 -> 0 (x = x0, z in [z0,z1], traveling toward -Z)
    {
        const double xc = x0;
        const double zc = clampd(z_m, z0, z1);
        const double denom = dz;
        double u = 0.0;
        if (std::abs(denom) > 1e-24) {
            // u=0 at z=z1, u=1 at z=z0
            u = (z1 - zc) / denom;
            u = clampd(u, 0.0, 1.0);
        }
        cands[3].s_m_raw = len0 + len1 + len2 + u * len3;
        cands[3].pos_room_m = make_v3(xc, geo_.y_m, zc);
        cands[3].d2 = sqr(x_m - xc) + sqr(z_m - zc);
        cands[3].seg = 3;
    }

    // Find best distance.
    int best = 0;
    double best_d2 = cands[0].d2;
    for (int i = 1; i < CeilingRailGeometry::kNumSegments; ++i) {
        if (cands[i].d2 < best_d2) {
            best_d2 = cands[i].d2;
            best = i;
        }
    }

    // Gather ties within tolerance.
    const double tol = 1e-12 * (1.0 + best_d2);

    // If multiple candidates tie, break ties deterministically using s_hint_m.
    // This is primarily relevant near corners.
    const double hint_w = wrapS(s_hint_m);

    double best_cost = std::numeric_limits<double>::infinity();
    int best_idx = best;

    for (int i = 0; i < CeilingRailGeometry::kNumSegments; ++i) {
        if (cands[i].d2 > best_d2 + tol) {
            continue;
        }

        // Candidate s for delta computation uses wrapped value.
        // Note: the closing corner may yield s==perimeter; wrapS maps it to 0.
        const double s_cand_w = wrapS(cands[i].s_m_raw);
        const double cost = std::abs(shortestDeltaS(hint_w, s_cand_w));

        if (cost < best_cost - 1e-15) {
            best_cost = cost;
            best_idx = i;
        } else if (std::abs(cost - best_cost) <= 1e-15) {
            // Stable fallback: prefer lower segment index.
            if (cands[i].seg < cands[best_idx].seg) {
                best_idx = i;
            }
        }
    }

    const Cand& c = cands[best_idx];

    out.s_m = wrapS(c.s_m_raw);
    out.pos_room_m = c.pos_room_m;
    out.dist_m = std::sqrt(std::max(0.0, c.d2));
    out.segment_idx = c.seg;

    return out;
}

} // namespace world
} // namespace vfep
