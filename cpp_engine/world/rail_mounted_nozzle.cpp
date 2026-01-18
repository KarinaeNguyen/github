#include "rail_mounted_nozzle.h"

#include <algorithm>
#include <cmath>

namespace vfep {
namespace world {

// Numeric stability thresholds
static constexpr double kEpsilon = 1e-12;
static constexpr double kMinVecLen = 1e-12;
static constexpr double kPI = 3.14159265358979323846;

static inline double clamp(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static inline double clamp01(double x) {
    return (x < 0.0) ? 0.0 : (x > 1.0) ? 1.0 : x;
}

static inline Vec3d v3(double x, double y, double z) { return Vec3d{x,y,z}; }

static inline Vec3d add(const Vec3d& a, const Vec3d& b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3d sub(const Vec3d& a, const Vec3d& b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline Vec3d mul(const Vec3d& a, double s)       { return v3(a.x*s,   a.y*s,   a.z*s); }

static inline double dot(const Vec3d& a, const Vec3d& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

static inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
    return v3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

static inline double len(const Vec3d& a) { return std::sqrt(dot(a,a)); }

static inline Vec3d norm(const Vec3d& a) {
    const double l = len(a);
    return (l > kMinVecLen) ? mul(a, 1.0 / l) : v3(0.0, 0.0, 0.0);
}

// Rodrigues rotation about unit axis.
static inline Vec3d rotate_axis_angle(const Vec3d& v, const Vec3d& axis_unit, double ang_rad) {
    const Vec3d a = norm(axis_unit);
    if (len(a) < kMinVecLen) return v;
    const double c = std::cos(ang_rad);
    const double s = std::sin(ang_rad);
    // v' = v*c + (a×v)*s + a*(a·v)*(1-c)
    return add(
        add(mul(v, c), mul(cross(a, v), s)),
        mul(a, dot(a, v) * (1.0 - c))
    );
}

// Map s in [0..1] to a point on the rail perimeter and a tangent direction.
// Uses CeilingRailGeometry corners order consistent with main_vis drawing:
// 0->1->2->3->0.
static bool rail_point_and_tangent(const CeilingRailGeometry& g, double s_0_1, Vec3d& out_p, Vec3d& out_tan_unit) {
    const Vec3d c0 = g.corners_room_m[0];
    const Vec3d c1 = g.corners_room_m[1];
    const Vec3d c2 = g.corners_room_m[2];
    const Vec3d c3 = g.corners_room_m[3];

    const Vec3d e0 = sub(c1, c0);
    const Vec3d e1 = sub(c2, c1);
    const Vec3d e2 = sub(c3, c2);
    const Vec3d e3 = sub(c0, c3);

    const double L0 = len(e0);
    const double L1 = len(e1);
    const double L2 = len(e2);
    const double L3 = len(e3);

    const double L = L0 + L1 + L2 + L3;
    if (!(L > 1e-9)) return false;

    double dist = clamp01(s_0_1) * L;

    auto seg = [&](const Vec3d& a, const Vec3d& d, double Li) -> bool {
        if (Li <= 1e-12) return false;
        const double t = clamp(dist / Li, 0.0, 1.0);
        out_p = add(a, mul(d, t));
        out_tan_unit = norm(d);
        return true;
    };

    if (dist <= L0) return seg(c0, e0, L0);
    dist -= L0;
    if (dist <= L1) return seg(c1, e1, L1);
    dist -= L1;
    if (dist <= L2) return seg(c2, e2, L2);
    dist -= L2;
    return seg(c3, e3, L3);
}

void RailMountedNozzle::recompute(const Inputs& in) {
    valid_ = false;

    if (in.ceiling_rail == nullptr) return;
    if (!in.ceiling_rail->isValid()) return;

    const auto& g = in.ceiling_rail->geometry();

    Vec3d rail_p, rail_tan;
    if (!rail_point_and_tangent(g, in.s_0_1, rail_p, rail_tan)) return;
    if (len(rail_tan) < 1e-12) rail_tan = v3(1.0, 0.0, 0.0);

    const Vec3d up = v3(0.0, 1.0, 0.0);

    // Nozzle position is rail position minus drop along +Y.
    const Vec3d nozzle_p = sub(rail_p, mul(up, cfg_.nozzle_drop_from_rail_m));

    // Base forward is tangent; build right-handed frame.
    Vec3d fwd0 = norm(rail_tan);
    Vec3d right0 = norm(cross(up, fwd0));
    if (len(right0) < 1e-12) right0 = v3(1.0, 0.0, 0.0);

    // Apply yaw about world up, then pitch about local right (after yaw).
    // Clamp angles to prevent numerical issues from extreme values
    double yaw_rad   = in.yaw_deg   * (kPI / 180.0);
    double pitch_rad = in.pitch_deg * (kPI / 180.0);
    
    // Normalize angles to [-pi, pi] range for stability
    const double k2PI = 2.0 * kPI;
    while (yaw_rad > kPI) yaw_rad -= k2PI;
    while (yaw_rad < -kPI) yaw_rad += k2PI;
    while (pitch_rad > kPI) pitch_rad -= k2PI;
    while (pitch_rad < -kPI) pitch_rad += k2PI;

    const Vec3d fwd1 = rotate_axis_angle(fwd0, up, yaw_rad);
    Vec3d right1 = norm(cross(up, fwd1));
    if (len(right1) < 1e-12) right1 = right0;

    const Vec3d fwd2 = rotate_axis_angle(fwd1, right1, pitch_rad);
    const Vec3d spray_dir = norm(fwd2);

    pose_.rail_pos_room_m        = rail_p;
    pose_.rail_tangent_unit_room = fwd0;

    pose_.nozzle_pos_room_m      = nozzle_p;
    pose_.spray_dir_unit_room    = (len(spray_dir) > kMinVecLen) ? spray_dir : fwd1;

    pose_.up_unit_room           = up;
    pose_.right_unit_room        = right1;
    pose_.forward_unit_room      = pose_.spray_dir_unit_room;

    valid_ = true;
}

} // namespace world
} // namespace vfep
