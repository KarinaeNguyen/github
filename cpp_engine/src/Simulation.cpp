// Simulation.cpp (patched: seedAmbient fail-safe + UI-friendly "safe-for-N-seconds" termination)
// Notes:
// - Implements safeHold_s_ (member in Simulation.h) for robust safe-for-N-seconds conclusion.
// - Recomputes inert_kgm3_ after ventilation so observations reflect end-of-step reactor state.

#include "Simulation.h"
#include "Aerodynamics.h"
#include "Constants.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace vfep {

namespace {

static std::uint32_t xorshift32(std::uint32_t& s) {
    // Deterministic, fast, portable.
    s ^= (s << 13);
    s ^= (s >> 17);
    s ^= (s << 5);
    return s;
}

static double u01(std::uint32_t& s) {
    // [0,1)
    return (double)xorshift32(s) / 4294967296.0; // 2^32
}

// --------------------
// Tunable constants
// --------------------
constexpr double kT_amb_K = 295.15;

constexpr double kAmbient_yO2  = 0.2095;
constexpr double kAmbient_yCO2 = 0.00042;
constexpr double kAmbient_yH2O = 0.0100;

constexpr double kPyrolysisStep_kgps = 0.01;

constexpr double kConclude_HRR_W   = 100.0;
constexpr double kConclude_T_C     = 35.0;
constexpr double kConclude_fuel_kg = 0.5;

// Additional UI-friendly conclude rule: if the system has been "safe" for long enough,
// conclude even if solid fuel remains (common when suppression succeeds early).
constexpr double kConcludeSafe_HRR_W   = 50.0;   // stricter than kConclude_HRR_W
constexpr double kConcludeSafe_T_C     = 40.0;   // modest safety threshold
constexpr double kConcludeSafeHold_s   = 60.0;   // must remain safe continuously

constexpr double kRewardSafe_T_C      = 100.0;
constexpr double kRewardSafeBonus     = 10.0;
constexpr double kRewardUnsafePenalty = -10.0;
constexpr double kRewardWasteCoeff    = 2.0;

// A tiny number for guarding divides / comparisons
constexpr double kEps = 1e-12;

static inline double clamp01(double x) {
    if (!std::isfinite(x)) return 0.0;
    return std::clamp(x, 0.0, 1.0);
}


static inline std::uint32_t fnv1a32_update(std::uint32_t h, const void* data, std::size_t len) {
    const std::uint8_t* p = reinterpret_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < len; ++i) {
        h ^= static_cast<std::uint32_t>(p[i]);
        h *= 16777619u;
    }
    return h;
}

static inline std::uint32_t fnv1a32_begin() { return 2166136261u; }

static inline std::uint32_t fnv1a32_add_u32(std::uint32_t h, std::uint32_t v) {
    return fnv1a32_update(h, &v, sizeof(v));
}
static inline std::uint32_t fnv1a32_add_i32(std::uint32_t h, std::int32_t v) {
    return fnv1a32_update(h, &v, sizeof(v));
}
static inline std::uint32_t fnv1a32_add_f64(std::uint32_t h, double v) {
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(v), "unexpected double size");
    std::memcpy(&bits, &v, sizeof(v));
    return fnv1a32_update(h, &bits, sizeof(bits));
}

static inline std::uint32_t crc32_update(std::uint32_t crc, const void* data, std::size_t len) {
    static bool table_init = false;
    static std::uint32_t table[256];
    if (!table_init) {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) {
                c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            }
            table[i] = c;
        }
        table_init = true;
    }
    const std::uint8_t* p = reinterpret_cast<const std::uint8_t*>(data);
    std::uint32_t c = crc ^ 0xFFFFFFFFu;
    for (std::size_t i = 0; i < len; ++i) {
        c = table[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    }
    return c ^ 0xFFFFFFFFu;
}

static inline std::uint32_t crc32_add_f32(std::uint32_t crc, float v) {
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(v), "unexpected float size");
    std::memcpy(&bits, &v, sizeof(v));
    return crc32_update(crc, &bits, sizeof(bits));
}

static inline double utilizationU(double exposure_kg, double k_util_1_per_kg) {
    if (!std::isfinite(exposure_kg) || exposure_kg <= 0.0) return 0.0;
    if (!std::isfinite(k_util_1_per_kg) || k_util_1_per_kg <= 0.0) return 0.0;
    const double u = 1.0 - std::exp(-k_util_1_per_kg * exposure_kg);
    return clamp01(u);
}

static inline double knockdownTargetFromEffectiveExposure(double effective_exposure_kg,
                                                         double EC50_adj_kg,
                                                         double hill,
                                                         double eps) {
    if (!std::isfinite(effective_exposure_kg) || effective_exposure_kg <= 0.0) return 0.0;
    if (!std::isfinite(EC50_adj_kg) || EC50_adj_kg <= eps) return 1.0;
    if (!std::isfinite(hill) || hill <= 0.0) hill = 1.0;

    const double ratio = EC50_adj_kg / (effective_exposure_kg + eps);
    const double p = std::pow(ratio, hill);
    const double kd = 1.0 / (1.0 + p);
    return clamp01(kd);
}


static inline double hillSaturation(double exposure_kg, double half_kg, double n) {
    if (!std::isfinite(exposure_kg) || exposure_kg <= 0.0) return 0.0;
    if (!std::isfinite(half_kg) || half_kg <= kEps) return 1.0;
    if (!std::isfinite(n) || n <= 0.0) n = 1.0;

    const double x = exposure_kg / half_kg;
    const double xn = std::pow(x, n);
    // xn / (1 + xn) is monotonic and saturating
    const double y = xn / (1.0 + xn);
    return clamp01(y);
}

static inline double firstOrderLag(double current, double target, double dt, double tau_s) {
    if (!std::isfinite(current)) current = 0.0;
    if (!std::isfinite(target)) target = 0.0;
    current = clamp01(current);
    target  = clamp01(target);

    if (!std::isfinite(dt) || dt <= 0.0) return current;
    if (!std::isfinite(tau_s) || tau_s <= kEps) return target;

    const double alpha = 1.0 - std::exp(-dt / tau_s);
    return clamp01(current + alpha * (target - current));
}

// Phase 2C: split delivered flow across 4 fixed sectors using spray direction.
// Sector directions are shallow offsets around the nominal +Z axis (nozzle frame).
static inline std::array<double, 4> sectorWeightsFromSprayDir(const Vec3d& spray_dir_unit) {
    // Build 4 sector direction vectors (sx, sy, +1) then compute weights from dot products.
    constexpr double kSide = 0.35; // small lateral offset to create quadrants
    const std::array<Vec3d, 4> dirs = {{
        Vec3d{-kSide, -kSide, 1.0},
        Vec3d{+kSide, -kSide, 1.0},
        Vec3d{-kSide, +kSide, 1.0},
        Vec3d{+kSide, +kSide, 1.0},
    }};

    std::array<double, 4> w{{0.25, 0.25, 0.25, 0.25}};
    double sum = 0.0;

    for (int i = 0; i < 4; ++i) {
        Vec3d d = dirs[i];
        const double mag = std::sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
        if (mag > kEps) { d.x /= mag; d.y /= mag; d.z /= mag; }

        const double dot = (spray_dir_unit.x * d.x + spray_dir_unit.y * d.y + spray_dir_unit.z * d.z);
        // Use a powered ReLU dot to emphasize alignment without randomness.
        const double wi = std::pow(std::max(0.0, dot), 4.0);
        w[i] = (std::isfinite(wi) ? wi : 0.0);
        sum += w[i];
    }

    if (sum <= kEps) {
        return {{0.25, 0.25, 0.25, 0.25}};
    }
    for (double& wi : w) wi /= sum;
    return w;
}




// --------------------
// Phase 3 geometry helpers (deterministic, no ray tracing)
// --------------------

static inline bool rayAabbIntersect(const Vec3d& ro, const Vec3d& rd_unit,
                                    const AABBd& box, double& tEnter, double& tExit) {
    // Slab method; rd_unit must be normalized and finite.
    tEnter = 0.0;
    tExit  = 1e30;

    const double rox = ro.x - box.c.x;
    const double roy = ro.y - box.c.y;
    const double roz = ro.z - box.c.z;

    const double hx = std::max(0.0, box.h.x);
    const double hy = std::max(0.0, box.h.y);
    const double hz = std::max(0.0, box.h.z);

    auto slab = [&](double roA, double rdA, double hA) -> bool {
        if (!std::isfinite(roA) || !std::isfinite(rdA) || !std::isfinite(hA)) return false;

        if (std::abs(rdA) < kEps) {
            // Ray parallel to slab: must be inside
            return (roA >= -hA && roA <= hA);
        }
        const double inv = 1.0 / rdA;
        double t0 = (-hA - roA) * inv;
        double t1 = ( +hA - roA) * inv;
        if (t0 > t1) std::swap(t0, t1);
        tEnter = std::max(tEnter, t0);
        tExit  = std::min(tExit,  t1);
        return (tEnter <= tExit);
    };

    if (!slab(rox, rd_unit.x, hx)) return false;
    if (!slab(roy, rd_unit.y, hy)) return false;
    if (!slab(roz, rd_unit.z, hz)) return false;

    return (tExit >= 0.0);
}

static inline Vec3d safeNorm(const Vec3d& v) {
    const double m2 = v.x*v.x + v.y*v.y + v.z*v.z;
    if (!std::isfinite(m2) || m2 <= kEps*kEps) return {0.0, 0.0, 0.0};
    const double inv = 1.0 / std::sqrt(m2);
    return {v.x*inv, v.y*inv, v.z*inv};
}

static inline double lineOfAttackCurved_0_1(const Vec3d& spray_dir_unit,
                                            const Vec3d& surface_normal_unit,
                                            double loa_min_0_1,
                                            double loa_power) {
    // Attack is best when spray direction is aligned with -normal (i.e., pointing into the face).
    const Vec3d inc = {-spray_dir_unit.x, -spray_dir_unit.y, -spray_dir_unit.z};
    const double c = inc.x*surface_normal_unit.x + inc.y*surface_normal_unit.y + inc.z*surface_normal_unit.z;
    const double d = std::clamp(c, 0.0, 1.0);
    const double p = (std::isfinite(loa_power) && loa_power > 0.0) ? loa_power : 1.0;
    const double mn = std::clamp(loa_min_0_1, 0.0, 1.0);
    const double v = mn + (1.0 - mn) * std::pow(d, p);
    return clamp01(v);
}

static Species MakeSpecies(const char* name, double molarMass_kg_per_mol, double cp, Phase ph) {
    return Species{name, molarMass_kg_per_mol, cp, ph};
}

static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}

} // namespace

std::vector<Species> Simulation::buildDefaultSpecies() {
    std::vector<Species> sp;
    sp.reserve(7);

    sp.push_back(MakeSpecies("N2",    0.0280134, 29.0, Phase::Gas));
    sp.push_back(MakeSpecies("O2",    0.0319980, 29.4, Phase::Gas));
    sp.push_back(MakeSpecies("CO2",   0.0440095, 37.1, Phase::Gas));
    sp.push_back(MakeSpecies("H2O",   0.0180153, 34.0, Phase::Gas));
    sp.push_back(MakeSpecies("FUEL",  0.0140000, 40.0, Phase::Gas));
    sp.push_back(MakeSpecies("INERT", 0.0440095, 37.1, Phase::Gas));
    sp.push_back(MakeSpecies("INHIB", 0.1000000,  0.0, Phase::Aerosol));


    return sp;
}

void Simulation::seedAmbient(Reactor& r) {
    auto& n = r.moles();
    std::fill(n.begin(), n.end(), 0.0);

    const double V = r.config().volume_m3;
    const double T = r.temperatureK();

    // Fail-safe: if config is invalid, still seed *something* so the reactor has a valid gas state.
    // This prevents totalGasMoles/Cp from becoming zero and avoids flatlined visualizer plots.
    if (!isFinitePositive(V) || !isFinitePositive(T)) {
        if (idx_.iN2 >= 0 && idx_.iN2 < static_cast<int>(n.size())) {
            // Choose an arbitrary small but positive fill (1 mol) if we can't compute ideal-gas totals.
            n[idx_.iN2] = 1.0;
        }
        return;
    }

    const double nTot = (P_atm * V) / (R_universal * T);
    if (!std::isfinite(nTot) || nTot <= 0.0) {
        if (idx_.iN2 >= 0 && idx_.iN2 < static_cast<int>(n.size())) {
            n[idx_.iN2] = 1.0;
        }
        return;
    }

    // Clamp declared ambient fractions into [0, 1]
    double yO2  = std::clamp(kAmbient_yO2,  0.0, 1.0);
    double yCO2 = std::clamp(kAmbient_yCO2, 0.0, 1.0);
    double yH2O = std::clamp(kAmbient_yH2O, 0.0, 1.0);

    // Renormalize if their sum exceeds 1.0
    const double sum = yO2 + yCO2 + yH2O;
    double yN2 = 0.0;

    if (sum > 1.0) {
        const double inv = 1.0 / sum;
        yO2  *= inv;
        yCO2 *= inv;
        yH2O *= inv;
        yN2 = 0.0;
    } else {
        yN2 = 1.0 - sum;
    }

    // Write moles (indices assumed valid from buildDefaultSpecies)
    n[idx_.iO2]  = nTot * yO2;
    n[idx_.iCO2] = nTot * yCO2;
    n[idx_.iH2O] = nTot * yH2O;
    n[idx_.iN2]  = nTot * yN2;
}

Simulation::Simulation()
: idx_{0,1,2,3,4,5,6},
  reactor_(
      buildDefaultSpecies(),
      idx_,
      []{
          CombustionModel cm;
          cm.C = 1.0; cm.H = 2.0; cm.O = 0.0;
          cm.A = 2.0e6;
          cm.Ea = 8.0e4;
          cm.orderFuel = 1.0;
          cm.orderO2 = 1.0;
          cm.heatRelease_J_per_molFuel = 6.0e5;
          return cm;
      }()
  ),
  vent_(idx_),
  supp_(idx_) {

    {
        ReactorConfig rc;
        rc.volume_m3  = 120.0;
        rc.area_m2    = 180.0;
        rc.T_amb_K    = kT_amb_K;
        rc.h_W_m2K    = 10.0;
        rc.emissivity = 0.85;
        reactor_.setConfig(rc);
    }



    {
        VentilationConfig vc;
        vc.ACH        = 3.0;
        vc.T_supply_K = kT_amb_K;
        vent_.setConfig(vc);
    }

    {
        SuppressionConfig sc;
        sc.enabled            = false;
        sc.mdot_total_kgps    = 0.30;
        sc.frac_inhibitor     = 0.25;
        sc.frac_inert         = 0.75;
        sc.cooling_W_per_kgps = 250000.0;
        supp_.setConfig(sc);
    }


// Phase 3B: default agent profile (science-ready placeholders; tunable in HUD).
setAgent(AgentType::CleanAgent);

    resetToDataCenterRackScenario();
}

void Simulation::resetToDataCenterRackScenario() {
    concluded_ = false;
    ignited_   = false;

    fuelSolid_kg_      = 50.0;
    pyrolysis_kgps_    = 0.0;
    pyrolysisMax_kgps_ = 0.06;

    reactor_.setTemperatureK(kT_amb_K);
    seedAmbient(reactor_);

    supp_.resetTank(200.0);
    {
        auto sc = supp_.config();
        sc.enabled = false;
        supp_.setConfig(sc);
    }

    // Li-ion runaway configuration
    LiIonConfig bc;
    bc.mass_pack_kg      = 8.0;
    bc.T_onset_K         = 423.15;
    bc.T_full_K          = 773.15;
    bc.E_total_J_per_kg  = 6.0e6;
    bc.tau_s             = 35.0;
    bc.ventGas_kg_per_kg = 0.12;

    liion_.configure(bc);
    liion_.setEnabled(true);
    liionHeat_W_    = 0.0;
    liionVent_kgps_ = 0.0;

    lastHRR_W_       = 0.0;
    inhib_kgm3_      = 0.0;
    inert_kgm3_      = 0.0;
    agent_mdot_kgps_ = 0.0;

    // Phase 2A truth telemetry defaults
    vfep_rpm_ = 0.0;
    hit_efficiency_0_1_ = 0.0;
    spray_dir_unit_ = {0.0, 0.0, 1.0};
    draft_vel_mps_  = {0.0, 0.0, 0.0};
    jet_momentum_N_ = 0.0;
    draft_drag_N_   = 0.0;

    // Phase 2B/2C suppression state defaults
    exposure_kg_ = 0.0;
    knockdown_target_0_1_ = 0.0;
    knockdown_0_1_ = 0.0;
    raw_HRR_W_ = 0.0;
    effective_HRR_W_ = 0.0;
    suppression_regime_ = SuppressionRegime::None;
    sector_exposure_kg_.fill(0.0);
    sector_knockdown_target_0_1_.fill(0.0);
    sector_knockdown_0_1_.fill(0.0);
    sector_delivered_mdot_kgps_.fill(0.0);
    sector_occlusion_0_1_.fill(1.0);
    sector_line_attack_0_1_.fill(0.0);
    sector_net_delivered_mdot_kgps_.fill(0.0);
    sector_utilization_U_0_1_.fill(0.0);
    sector_effective_exposure_kg_.fill(0.0);
    sector_EC50_adj_kg_.fill(0.0);
    utilization_U_0_1_ = 0.0;
    effective_exposure_kg_ = 0.0;
    EC50_adj_kg_ = 0.0;
    // NOTE: Phase 3B.1 signatures/telemetry are accumulated in a separate harness
    // sampler (do not reset here).

    // Phase 3B.1 harness state
    run_signatures_ = {};
    expected_signatures_ = {};
    latest_events_bits_ = 0;
    telemetry_head_ = 0;
    telemetry_count_ = 0;
    telemetry_next_t_s_ = 0.0;
    prev_occluded_any_ = false;
    prev_kd_ge_0p5_ = false;
    prev_kd_ge_0p9_ = false;
    prev_hrr_below_100kW_ = false;
    prev_exposure_kg_ = 0.0;
    prev_effective_exposure_kg_ = 0.0;
    prev_kd_target_0_1_ = 0.0;

    // Phase 3A stabilization + audit telemetry state defaults
    sector_raw_delivered_mdot_kgps_.fill(0.0);
    sector_shield_0_1_.fill(1.0);
    occ_stable_0_1_.fill(1.0);
    shield_stable_0_1_.fill(1.0);
    loa_smooth_0_1_.fill(0.0);
    occ_enter_count_.fill(0);
    occ_exit_count_.fill(0);
    shield_enter_count_.fill(0);
    shield_exit_count_.fill(0);
    fully_blocked_hold_s_ = 0.0;
    glancing_hold_s_ = 0.0;

    scenario_time_s_ = 0.0;
    nozzle_sweep_enabled_ = false;
    num_obstacles_ = 0;

    // Match visualizer nozzle axis (Phase 1) but keep it sim-owned for determinism.
    // (Will be normalized implicitly by aero module.)
    nozzle_dir_unit_scenario_ = {0.7, -0.15, 0.7};
    nozzle_pos_m_ = {-2.0, 1.5, -2.0};
    hotspot_pos_m_ = {0.0, 0.6, 0.41};
    ignition_seeded_ = false;
    ignition_seed_u32_ = 0u;
    nozzle_dir_unit_ = nozzle_dir_unit_scenario_;

    safeHold_s_ = 0.0;
}

void Simulation::setAgent(AgentType a) {
    agent_type_ = a;

    // Placeholder calibration table (explicit units; tune as needed).
    // NOTE: These are not claims; they are starting points for calibration.
    switch (a) {
    case AgentType::CleanAgent:
        agent_profile_.EC50_kg = 0.80;
        agent_profile_.hill = 1.6;
        agent_profile_.k_util_1_per_kg = 1.2;
        agent_profile_.potency = 1.0;
        break;
    case AgentType::DryChemical:
        agent_profile_.EC50_kg = 0.25;
        agent_profile_.hill = 2.4;
        agent_profile_.k_util_1_per_kg = 4.0;
        agent_profile_.potency = 1.0;
        break;
    case AgentType::CO2:
    default:
        agent_profile_.EC50_kg = 1.10;
        agent_profile_.hill = 1.3;
        agent_profile_.k_util_1_per_kg = 0.8;
        agent_profile_.potency = 1.0;
        break;
    }

    // Keep signatures deterministic and audit-ready.
    std::uint32_t h = fnv1a32_begin();
    h = fnv1a32_add_i32(h, static_cast<std::int32_t>(agent_type_));
    h = fnv1a32_add_f64(h, agent_profile_.EC50_kg);
    h = fnv1a32_add_f64(h, agent_profile_.hill);
    h = fnv1a32_add_f64(h, agent_profile_.k_util_1_per_kg);
    h = fnv1a32_add_f64(h, agent_profile_.potency);
    h = fnv1a32_add_f64(h, scenario_factors_.temp_factor);
    h = fnv1a32_add_f64(h, scenario_factors_.vent_factor);
    h = fnv1a32_add_f64(h, scenario_factors_.fuel_factor);
    run_param_hash_u32_ = h;
}

void Simulation::resetToScenario(DemoScenario s, AgentType a) {
    setAgent(a);
    resetToScenario(s);
}

void Simulation::enableCalibrationMode(bool enabled) {
    calibration_mode_ = enabled;
    calibration_elapsed_s_ = 0.0;
    telemetry_crc_u32_ = 0;
    // Phase 3B.1 harness signatures (shared across calibration/verification)
    run_signatures_ = {};
    expected_signatures_ = {};
    latest_events_bits_ = 0;
    telemetry_head_ = 0;
    telemetry_count_ = 0;
    telemetry_next_t_s_ = 0.0;
    telemetry_v2_head_ = 0;
    telemetry_v2_count_ = 0;
    last_profile_ = {};

    if (calibration_mode_) {
        // Deterministic calibration setup:
        // - fixed geometry (no obstacles)
        // - fixed nozzle pose (no sweep)
        // - fixed mass delivery (via mdot + deterministic aero)
        resetToScenario(DemoScenario::DirectVsGlance, agent_type_);
        num_obstacles_ = 0;
        nozzle_sweep_enabled_ = false;
        // Aim directly into the rack proxy (stable LoA).
        nozzle_pos_m_ = {-2.0, 1.5, -2.0};
        nozzle_dir_unit_scenario_ = {0.0, 0.0, 1.0};
        // Use a stable, explicit discharge for calibration harness runs.
        agent_mdot_kgps_ = 0.30;
        vfep_rpm_ = 0.0;
        // Scenario factors neutral by default in calibration harness; user can override in code/HUD later.
        scenario_factors_.temp_factor = 1.0;
        scenario_factors_.vent_factor = 1.0;
        scenario_factors_.fuel_factor = 1.0;

        // Recompute signature now that all deterministic calibration parameters are set.
        std::uint32_t h = fnv1a32_begin();
        h = fnv1a32_add_i32(h, static_cast<std::int32_t>(agent_type_));
        h = fnv1a32_add_f64(h, agent_profile_.EC50_kg);
        h = fnv1a32_add_f64(h, agent_profile_.hill);
        h = fnv1a32_add_f64(h, agent_profile_.k_util_1_per_kg);
        h = fnv1a32_add_f64(h, agent_profile_.potency);
        h = fnv1a32_add_f64(h, scenario_factors_.temp_factor);
        h = fnv1a32_add_f64(h, scenario_factors_.vent_factor);
        h = fnv1a32_add_f64(h, scenario_factors_.fuel_factor);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.x);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.y);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.z);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.x);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.y);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.z);
        h = fnv1a32_add_f64(h, agent_mdot_kgps_);
        h = fnv1a32_add_f64(h, vent_.config().ACH);
        run_param_hash_u32_ = h;
        run_signatures_.run_param_hash_u32 = run_param_hash_u32_;
    } else {
        // Leaving calibration mode keeps current scenario state; signatures remain computed from current parameters.
        setAgent(agent_type_);
    }
}

void Simulation::enableVerificationMode(bool enabled) {
    verification_mode_ = enabled;
    // Reset harness outputs (do not mutate core state)
    run_signatures_ = {};
    expected_signatures_ = {};
    latest_events_bits_ = 0;
    telemetry_head_ = 0;
    telemetry_count_ = 0;
    telemetry_next_t_s_ = 0.0;
    telemetry_v2_head_ = 0;
    telemetry_v2_count_ = 0;
    last_profile_ = {};
    prev_occluded_any_ = false;
    prev_kd_ge_0p5_ = false;
    prev_kd_ge_0p9_ = false;
    prev_hrr_below_100kW_ = false;
    prev_exposure_kg_ = 0.0;
    prev_effective_exposure_kg_ = 0.0;
    prev_kd_target_0_1_ = 0.0;
}

void Simulation::enableTelemetryV2(bool enabled) {
    telemetry_v2_enabled_ = enabled;
    telemetry_v2_head_ = 0;
    telemetry_v2_count_ = 0;
}

int Simulation::getTelemetrySamplesV2(TelemetrySampleV2* out_ptr, int cap) const {
    if (!out_ptr || cap <= 0) return 0;
    const int n = std::min<int>(telemetry_v2_count_, cap);
    int idx = (telemetry_v2_head_ - telemetry_v2_count_);
    while (idx < 0) idx += kTelemetryV2Capacity_;
    for (int i = 0; i < n; ++i) {
        out_ptr[i] = telemetry_v2_rb_[(idx + i) % kTelemetryV2Capacity_];
    }
    return n;
}

void Simulation::setControlInputs(const ControlInputsV1& in) {
    // Deterministic contract enforcement: accept only matching version/size.
    if (in.version_u32 != 1u || in.size_bytes_u32 != sizeof(ControlInputsV1)) {
        return;
    }
    ControlInputsV1 tmp = in;
    // Clamp to safe ranges; no side effects on physics if defaults are used.
    if (!std::isfinite(tmp.mdot_scale) || tmp.mdot_scale < 0.0) tmp.mdot_scale = 0.0;
    if (tmp.mdot_scale > 1000.0) tmp.mdot_scale = 1000.0;
    // Mask only legacy sectors for Phase 3CA baseline.
    tmp.sector_enable_mask_u32 &= 0xFu;
    control_inputs_ = tmp;
}

void Simulation::enableProfiling(bool enabled) {
    profiling_enabled_ = (enabled && !verification_mode_);
    last_profile_ = {};
}

bool Simulation::getLastProfileSample(ProfileSampleV1* out) const {
    if (!out) return false;
    *out = last_profile_;
    return true;
}

RunSignatures Simulation::getRunSignatures() const { return run_signatures_; }

RunSignatures Simulation::getLastExpectedSignatures() const { return expected_signatures_; }

std::uint32_t Simulation::getLatestEvents() {
    const std::uint32_t out = latest_events_bits_;
    latest_events_bits_ = 0;
    return out;
}

int Simulation::getTelemetrySamples(TelemetrySampleV1* out_ptr, int cap) const {
    if (!out_ptr || cap <= 0) return 0;
    const int n = std::min<int>(telemetry_count_, cap);
    // Oldest sample index = head - count (mod capacity)
    int idx = (telemetry_head_ - telemetry_count_);
    while (idx < 0) idx += kTelemetryCapacity_;
    for (int i = 0; i < n; ++i) {
        out_ptr[i] = telemetry_rb_[(idx + i) % kTelemetryCapacity_];
    }
    return n;
}

static inline std::uint32_t fnv_hash_text_u32(const char* s) {
    if (!s) return 0;
    std::uint32_t h = fnv1a32_begin();
    h = fnv1a32_update(h, s, std::strlen(s));
    return h;
}

int Simulation::exportConfigText(char* buf, int cap) const {
    if (!buf || cap <= 0) return 0;

    // Build versioned contracts with explicit units.
    SimConfigV1 sc;
    sc.fixed_dt_s = 0.0; // runtime-selected by harness (verification) or UI
    sc.telemetry_dt_s = telemetry_dt_s_;
    sc.exposure_half_kg = exposure_half_kg_;
    sc.exposure_hill_n = exposure_hill_n_;
    sc.tau_knockdown_rise_s = tau_knockdown_rise_s_;
    sc.tau_knockdown_fall_s = tau_knockdown_fall_s_;
    sc.enable_recovery_u32 = enable_recovery_ ? 1u : 0u;
    sc.tau_exposure_decay_s = tau_exposure_decay_s_;
    sc.fnv_hash_u32 = 0;

    AgentConfigV1 ac;
    ac.agent_type_i32 = static_cast<std::int32_t>(agent_type_);
    ac.EC50_kg = agent_profile_.EC50_kg;
    ac.hill = agent_profile_.hill;
    ac.k_util_1_per_kg = agent_profile_.k_util_1_per_kg;
    ac.potency = agent_profile_.potency;
    ac.fnv_hash_u32 = 0;

    ScenarioConfigV1 cc;
    cc.demo_scenario_i32 = static_cast<std::int32_t>(scenario_);
    cc.nozzle_pos_m = nozzle_pos_m_;
    cc.nozzle_dir_unit = nozzle_dir_unit_scenario_;
    cc.nozzle_sweep_enabled_u32 = nozzle_sweep_enabled_ ? 1u : 0u;
    cc.sweep_amp_deg = sweep_amp_deg_;
    cc.sweep_freq_hz = sweep_freq_hz_;
    cc.mdot_cmd_kgps = agent_mdot_kgps_;
    cc.ACH_1_per_h = vent_.config().ACH;
    cc.temp_factor = scenario_factors_.temp_factor;
    cc.vent_factor = scenario_factors_.vent_factor;
    cc.fuel_factor = scenario_factors_.fuel_factor;
    cc.fnv_hash_u32 = 0;

    // Deterministic field hashing (explicit list, fixed order)
    auto hash_sim = [&]() {
        std::uint32_t h = fnv1a32_begin();
        h = fnv1a32_add_u32(h, sc.version_u32);
        h = fnv1a32_add_u32(h, sc.size_bytes_u32);
        h = fnv1a32_add_f64(h, sc.fixed_dt_s);
        h = fnv1a32_add_f64(h, sc.telemetry_dt_s);
        h = fnv1a32_add_f64(h, sc.exposure_half_kg);
        h = fnv1a32_add_f64(h, sc.exposure_hill_n);
        h = fnv1a32_add_f64(h, sc.tau_knockdown_rise_s);
        h = fnv1a32_add_f64(h, sc.tau_knockdown_fall_s);
        h = fnv1a32_add_u32(h, sc.enable_recovery_u32);
        h = fnv1a32_add_f64(h, sc.tau_exposure_decay_s);
        return h;
    };
    auto hash_agent = [&]() {
        std::uint32_t h = fnv1a32_begin();
        h = fnv1a32_add_u32(h, ac.version_u32);
        h = fnv1a32_add_u32(h, ac.size_bytes_u32);
        h = fnv1a32_add_i32(h, ac.agent_type_i32);
        h = fnv1a32_add_f64(h, ac.EC50_kg);
        h = fnv1a32_add_f64(h, ac.hill);
        h = fnv1a32_add_f64(h, ac.k_util_1_per_kg);
        h = fnv1a32_add_f64(h, ac.potency);
        return h;
    };
    auto hash_scn = [&]() {
        std::uint32_t h = fnv1a32_begin();
        h = fnv1a32_add_u32(h, cc.version_u32);
        h = fnv1a32_add_u32(h, cc.size_bytes_u32);
        h = fnv1a32_add_i32(h, cc.demo_scenario_i32);
        h = fnv1a32_add_f64(h, cc.nozzle_pos_m.x);
        h = fnv1a32_add_f64(h, cc.nozzle_pos_m.y);
        h = fnv1a32_add_f64(h, cc.nozzle_pos_m.z);
        h = fnv1a32_add_f64(h, cc.nozzle_dir_unit.x);
        h = fnv1a32_add_f64(h, cc.nozzle_dir_unit.y);
        h = fnv1a32_add_f64(h, cc.nozzle_dir_unit.z);
        h = fnv1a32_add_u32(h, cc.nozzle_sweep_enabled_u32);
        h = fnv1a32_add_f64(h, cc.sweep_amp_deg);
        h = fnv1a32_add_f64(h, cc.sweep_freq_hz);
        h = fnv1a32_add_f64(h, cc.mdot_cmd_kgps);
        h = fnv1a32_add_f64(h, cc.ACH_1_per_h);
        h = fnv1a32_add_f64(h, cc.temp_factor);
        h = fnv1a32_add_f64(h, cc.vent_factor);
        h = fnv1a32_add_f64(h, cc.fuel_factor);
        return h;
    };

    sc.fnv_hash_u32 = hash_sim();
    ac.fnv_hash_u32 = hash_agent();
    cc.fnv_hash_u32 = hash_scn();

    // Stable export text (deterministic order)
    int n = 0;
    auto app = [&](const char* fmt, auto... args) {
        if (n >= cap) return;
        const int w = std::snprintf(buf + n, (size_t)(cap - n), fmt, args...);
        if (w > 0) n += std::min(w, cap - n);
    };

    app("SimConfigV1\n");
    app("  version_u32=%u\n", sc.version_u32);
    app("  size_bytes_u32=%u\n", sc.size_bytes_u32);
    app("  fixed_dt_s=%.6f\n", sc.fixed_dt_s);
    app("  telemetry_dt_s=%.6f\n", sc.telemetry_dt_s);
    app("  exposure_half_kg=%.9g\n", sc.exposure_half_kg);
    app("  exposure_hill_n=%.9g\n", sc.exposure_hill_n);
    app("  tau_knockdown_rise_s=%.9g\n", sc.tau_knockdown_rise_s);
    app("  tau_knockdown_fall_s=%.9g\n", sc.tau_knockdown_fall_s);
    app("  enable_recovery_u32=%u\n", sc.enable_recovery_u32);
    app("  tau_exposure_decay_s=%.9g\n", sc.tau_exposure_decay_s);
    app("  fnv_hash_u32=0x%08X\n", sc.fnv_hash_u32);

    app("AgentConfigV1\n");
    app("  version_u32=%u\n", ac.version_u32);
    app("  size_bytes_u32=%u\n", ac.size_bytes_u32);
    app("  agent_type_i32=%d\n", ac.agent_type_i32);
    app("  EC50_kg=%.9g\n", ac.EC50_kg);
    app("  hill=%.9g\n", ac.hill);
    app("  k_util_1_per_kg=%.9g\n", ac.k_util_1_per_kg);
    app("  potency=%.9g\n", ac.potency);
    app("  fnv_hash_u32=0x%08X\n", ac.fnv_hash_u32);

    app("ScenarioConfigV1\n");
    app("  version_u32=%u\n", cc.version_u32);
    app("  size_bytes_u32=%u\n", cc.size_bytes_u32);
    app("  demo_scenario_i32=%d\n", cc.demo_scenario_i32);
    app("  nozzle_pos_m=(%.6f,%.6f,%.6f)\n", cc.nozzle_pos_m.x, cc.nozzle_pos_m.y, cc.nozzle_pos_m.z);
    app("  nozzle_dir_unit=(%.6f,%.6f,%.6f)\n", cc.nozzle_dir_unit.x, cc.nozzle_dir_unit.y, cc.nozzle_dir_unit.z);
    app("  nozzle_sweep_enabled_u32=%u\n", cc.nozzle_sweep_enabled_u32);
    app("  sweep_amp_deg=%.6f\n", cc.sweep_amp_deg);
    app("  sweep_freq_hz=%.6f\n", cc.sweep_freq_hz);
    app("  mdot_cmd_kgps=%.6f\n", cc.mdot_cmd_kgps);
    app("  ACH_1_per_h=%.6f\n", cc.ACH_1_per_h);
    app("  temp_factor=%.6f\n", cc.temp_factor);
    app("  vent_factor=%.6f\n", cc.vent_factor);
    app("  fuel_factor=%.6f\n", cc.fuel_factor);
    app("  fnv_hash_u32=0x%08X\n", cc.fnv_hash_u32);

    // Convenience: full export hash for copy/paste audits.
    const std::uint32_t export_hash = fnv_hash_text_u32(buf);
    app("ExportTextHash(FNV-1a32)=0x%08X\n", export_hash);

    return std::min(n, cap);
}

bool Simulation::runVerificationTest(VerificationTestId id) {
    // Embedded golden vectors (no file I/O)
    struct Golden { std::uint32_t p, c, d; };
    static const Golden kGolden[] = {
        /*V0*/ {0xB36C622Cu, 0xDD583539u, 0x97368E72u},
        /*V1*/ {0x4C6CC74Du, 0xDD583539u, 0xC714D0B6u},
        /*V2*/ {0xD5C9F352u, 0x8E0DDB47u, 0x9AACDBF8u},
    };

    enableVerificationMode(true);
    calibration_mode_ = false;

    // Deterministic vector setup
    switch (id) {
        case VerificationTestId::V0:
            resetToScenario(DemoScenario::DirectVsGlance, AgentType::CleanAgent);
            num_obstacles_ = 0;
            nozzle_sweep_enabled_ = false;
            nozzle_pos_m_ = {-2.0, 1.5, -2.0};
            nozzle_dir_unit_scenario_ = {0.0, 0.0, 1.0};
            agent_mdot_kgps_ = 0.30;
            scenario_factors_.temp_factor = 1.0;
            scenario_factors_.vent_factor = 1.0;
            scenario_factors_.fuel_factor = 1.0;
            break;
        case VerificationTestId::V1:
            resetToScenario(DemoScenario::OcclusionWall, AgentType::CleanAgent);
            nozzle_sweep_enabled_ = false;
            nozzle_pos_m_ = {-2.0, 1.5, -2.0};
            nozzle_dir_unit_scenario_ = {0.0, 0.0, 1.0};
            agent_mdot_kgps_ = 0.30;
            break;
        case VerificationTestId::V2:
        default:
            resetToScenario(DemoScenario::Mixed, AgentType::CleanAgent);
            nozzle_pos_m_ = {-2.0, 1.5, -2.0};
            nozzle_dir_unit_scenario_ = {0.0, 0.0, 1.0};
            nozzle_sweep_enabled_ = true;
            sweep_amp_deg_ = 14.0;
            sweep_freq_hz_ = 0.50;
            agent_mdot_kgps_ = 0.30;
            break;
    }

    // Fixed cadence + fixed-step sim length
    telemetry_dt_s_ = 0.10;
    telemetry_next_t_s_ = 0.0;
    run_signatures_ = {};
    latest_events_bits_ = 0;
    telemetry_head_ = 0;
    telemetry_count_ = 0;

    // Parameter hash over effective parameters (explicit list, fixed order)
    {
        std::uint32_t h = fnv1a32_begin();
        h = fnv1a32_add_i32(h, static_cast<std::int32_t>(id));
        h = fnv1a32_add_i32(h, static_cast<std::int32_t>(agent_type_));
        h = fnv1a32_add_f64(h, agent_profile_.EC50_kg);
        h = fnv1a32_add_f64(h, agent_profile_.hill);
        h = fnv1a32_add_f64(h, agent_profile_.k_util_1_per_kg);
        h = fnv1a32_add_f64(h, agent_profile_.potency);
        h = fnv1a32_add_f64(h, scenario_factors_.temp_factor);
        h = fnv1a32_add_f64(h, scenario_factors_.vent_factor);
        h = fnv1a32_add_f64(h, scenario_factors_.fuel_factor);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.x);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.y);
        h = fnv1a32_add_f64(h, nozzle_pos_m_.z);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.x);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.y);
        h = fnv1a32_add_f64(h, nozzle_dir_unit_scenario_.z);
        h = fnv1a32_add_u32(h, nozzle_sweep_enabled_ ? 1u : 0u);
        h = fnv1a32_add_f64(h, sweep_amp_deg_);
        h = fnv1a32_add_f64(h, sweep_freq_hz_);
        h = fnv1a32_add_f64(h, agent_mdot_kgps_);
        h = fnv1a32_add_f64(h, vent_.config().ACH);
        run_signatures_.run_param_hash_u32 = h;
    }

    const double fixed_dt_s = 0.01;
    const double T_s = (id == VerificationTestId::V2) ? 12.0 : 8.0;
    const int steps = (int)std::ceil(T_s / fixed_dt_s);

    ignited_ = true;
    concluded_ = false;

    for (int i = 0; i < steps; ++i) {
        step(fixed_dt_s);
    }

    const std::size_t idx = static_cast<std::size_t>(id);
    expected_signatures_.run_param_hash_u32 = kGolden[idx].p;
    expected_signatures_.telemetry_crc_u32  = kGolden[idx].c;
    expected_signatures_.state_digest_u32   = kGolden[idx].d;

    const bool pass = (expected_signatures_.run_param_hash_u32 == run_signatures_.run_param_hash_u32)
        && (expected_signatures_.telemetry_crc_u32 == run_signatures_.telemetry_crc_u32)
        && (expected_signatures_.state_digest_u32 == run_signatures_.state_digest_u32);

    return pass;
}

void Simulation::resetToScenario(DemoScenario s) {
    resetToDataCenterRackScenario();
scenario_ = s;

// Phase 3B: deterministic, static scenario factors (pure scalars; no dynamics).
scenario_factors_ = {};
scenario_factors_.temp_factor = 1.0;
scenario_factors_.vent_factor = 1.0;
scenario_factors_.fuel_factor = 1.0;
// Minimal, explainable defaults per scenario (tunable, not claims).
switch (s) {
case DemoScenario::OcclusionWall:
    scenario_factors_.vent_factor = 1.10;
    break;
case DemoScenario::ShieldingStack:
    scenario_factors_.fuel_factor = 1.15;
    break;
case DemoScenario::Mixed:
    scenario_factors_.temp_factor = 1.05;
    scenario_factors_.vent_factor = 1.10;
    scenario_factors_.fuel_factor = 1.10;
    break;
case DemoScenario::DirectVsGlance:
default:
    break;
}

// Update parameter signature when scenario context changes.
setAgent(agent_type_);


    // Deterministic, investor-demo-safe scenario setups.
    // Keep these minimal and explainable: only nozzle pose and up to two AABBs.
    num_obstacles_ = 0;
    nozzle_sweep_enabled_ = false;
    sweep_freq_hz_ = 0.25;
    sweep_amp_deg_ = 12.0;

    // Shared anchors (mirror visualizer layout).
    nozzle_pos_m_ = {-2.0, 1.5, -2.0};

    switch (scenario_) {
        case DemoScenario::DirectVsGlance: {
            // No occluders. Operator can see direct vs glancing purely from LoA.
            nozzle_dir_unit_scenario_ = {0.55, -0.10, 0.83};
            break;
        }
        case DemoScenario::OcclusionWall: {
            // A wall that blocks the left side under small motion.
            nozzle_dir_unit_scenario_ = {0.65, -0.15, 0.74};
            obstacles_[0] = AABBd{{-0.75, 1.20, -0.35}, {0.25, 0.35, 0.90}};
            num_obstacles_ = 1;
            break;
        }
        case DemoScenario::ShieldingStack: {
            // No occluders; rely on proxy sector shielding (near sectors shadow far).
            nozzle_dir_unit_scenario_ = {0.40, -0.10, 0.91};
            break;
        }
        case DemoScenario::Mixed: {
            // One occluder plus deterministic sweep to exercise hysteresis.
            nozzle_dir_unit_scenario_ = {0.60, -0.10, 0.79};
            obstacles_[0] = AABBd{{-0.60, 1.10, 0.20}, {0.22, 0.40, 0.55}};
            num_obstacles_ = 1;
            nozzle_sweep_enabled_ = true;
            sweep_freq_hz_ = 0.20;
            sweep_amp_deg_ = 10.0;
            break;
        }
        default: {
            nozzle_dir_unit_scenario_ = {0.7, -0.15, 0.7};
            break;
        }
    }

    // Apply scenario nozzle pose to aero module input.
    nozzle_dir_unit_ = nozzle_dir_unit_scenario_;
    scenario_time_s_ = 0.0;
}

void Simulation::setNozzlePose(const Vec3d& pos_m, const Vec3d& dir_unit) {
    // Sim-owned; no side effects outside deterministic state.
    if (std::isfinite(pos_m.x) && std::isfinite(pos_m.y) && std::isfinite(pos_m.z)) {
        nozzle_pos_m_ = pos_m;
    }
    const Vec3d d = safeNorm(dir_unit);
    if (std::abs(d.x) > kEps || std::abs(d.y) > kEps || std::abs(d.z) > kEps) {
        nozzle_dir_unit_scenario_ = d;
        nozzle_dir_unit_ = d;
    }
}

void Simulation::setNozzleSweepEnabled(bool enabled) {
    nozzle_sweep_enabled_ = enabled;
}

void Simulation::commandIgniteOrIncreasePyrolysis() {
    const bool first_ignite = !ignited_;
    ignited_ = true;

    // On first ignition, choose a deterministic "random" hotspot near the rack.
    if (first_ignite) {
        if (!ignition_seeded_) {
            // Prefer deterministic run hash if available; otherwise fixed constant.
            ignition_seed_u32_ = (run_param_hash_u32_ != 0u) ? run_param_hash_u32_ : 0xC0FFEE01u;
            ignition_seeded_ = true;
        }

        std::uint32_t s = ignition_seed_u32_;

        // Rack-front ignition volume (meters, world frame).
        // Keep fire on the rack front face (z ~= +0.4 in the default scene).
        const double x = -0.50 + 1.00 * u01(s);  // [-0.50, +0.50]
        const double y =  0.40 + 1.20 * u01(s);  // [0.40, 1.60]
        const double z =  0.41;                  // slightly in front of rack face

        hotspot_pos_m_ = {x, y, z};

        // Persist advanced seed (stable for the run).
        ignition_seed_u32_ = s;
    }

    // Ensure pyrolysis stays within bounds and does not go negative.
    pyrolysis_kgps_ = std::clamp(pyrolysis_kgps_ + kPyrolysisStep_kgps, 0.0, pyrolysisMax_kgps_);
}

void Simulation::commandStartSuppression() {
    auto sc = supp_.config();
    sc.enabled = true;
    supp_.setConfig(sc);
}

void Simulation::step(double dt) {
    if (concluded_ && !verification_mode_) return;

    // Reject non-physical time steps.
    if (!isFinitePositive(dt)) {
        return;
    }

    // --------------------
    // Phase 3CA.1: deterministic profiling counters (observational only)
    // --------------------
    // - No wall-clock / OS time usage.
    // - Hard-disabled in VerificationMode.
    const bool profile_active = profiling_enabled_ && !verification_mode_;
    if (profile_active) {
        last_profile_ = {};
        last_profile_.version_u32 = 2u;
        last_profile_.size_bytes_u32 = sizeof(ProfileSampleV1);
        last_profile_.t_s = static_cast<float>(scenario_time_s_);
    }

    auto prof_add = [&](ProfileStage st, std::uint64_t work) {
        if (!profile_active) return;
        const std::size_t idx = static_cast<std::size_t>(st);
        if (idx < last_profile_.stage_work.size()) {
            last_profile_.stage_work[idx] += work;
        }
    };

    // --------------------
    // Phase 3A scenario timebase + optional deterministic nozzle sweep
    // --------------------
    scenario_time_s_ += dt;
    if (!std::isfinite(scenario_time_s_) || scenario_time_s_ < 0.0) scenario_time_s_ = 0.0;

    // Apply scenario nozzle pose to aero module input, optionally with a deterministic yaw sweep.
    {
        Vec3d d = nozzle_dir_unit_scenario_;
        if (nozzle_sweep_enabled_ && sweep_amp_deg_ > 0.0 && sweep_freq_hz_ > 0.0) {
            const double a = (sweep_amp_deg_ * 3.14159265358979323846 / 180.0)
                           * std::sin(2.0 * 3.14159265358979323846 * sweep_freq_hz_ * scenario_time_s_);
            const double ca = std::cos(a);
            const double sa = std::sin(a);
            // Yaw about +Y: x' = ca*x + sa*z, z' = -sa*x + ca*z
            d = {ca * d.x + sa * d.z, d.y, -sa * d.x + ca * d.z};
        }
        nozzle_dir_unit_ = safeNorm(d);
    }
    prof_add(ProfileStage::ScenarioSweep, 1);
    prof_add(ProfileStage::ScenarioSweep, 1);

    // --------------------
    // Pyrolysis -> FUEL gas
    // --------------------
    // Phase 2B/2C causal coupling: knockdown reduces volatilization / fuel feed deterministically.
    if (ignited_ && fuelSolid_kg_ > kEps && pyrolysis_kgps_ > 0.0) {
        const double kd = clamp01(knockdown_0_1_);
        const double m_raw = std::min(fuelSolid_kg_, pyrolysis_kgps_ * dt);

        // Suppression reduces effective pyrolysis feed smoothly (via kd smoothing).
        const double m_eff = std::clamp(m_raw * (1.0 - kd), 0.0, m_raw);

        fuelSolid_kg_ -= m_eff;

        const auto& sp = reactor_.species();
        const double Mfuel = sp[idx_.iFUEL].molarMass_kg_per_mol;

        if (Mfuel > kEps) {
            reactor_.addMoles(idx_.iFUEL, m_eff / Mfuel);
        }
    }
    prof_add(ProfileStage::PyrolysisFuel, 1);
    prof_add(ProfileStage::PyrolysisFuel, 1);

    // --------------------
    // Apply suppression
    // --------------------
    const double cooling_W =
        supp_.apply(dt, reactor_, inhib_kgm3_, inert_kgm3_, agent_mdot_kgps_);
    prof_add(ProfileStage::SuppressionApply, 1);
    prof_add(ProfileStage::SuppressionApply, 1);

    // Actuator truth telemetry
    vfep_rpm_ = supp_.vfepRPM();

    // --------------------
    // Phase 2A aero telemetry (truth) + Phase 2B/2C causal suppression
    // --------------------
    // Deterministic cross-draft model: magnitude scales with ACH, direction fixed along +X.
    const double ach = vent_.config().ACH;
    const double draft_mag_mps = std::max(0.0, 0.35 * ach);
    draft_vel_mps_ = {draft_mag_mps, 0.0, 0.0};

    // Compute truth spray direction + hit efficiency.
    {
        AeroInputs ai;
        ai.mdot_kgps = agent_mdot_kgps_;
        ai.vfep_rpm = vfep_rpm_;
        ai.draft_vel_mps = draft_vel_mps_;
        ai.nozzle_dir_unit = nozzle_dir_unit_;

        const AeroOutputs ao = computeAerodynamics(ai);
        hit_efficiency_0_1_ = clamp01(ao.hit_efficiency_0_1);
        spray_dir_unit_ = ao.spray_dir_unit;
        jet_momentum_N_ = ao.jet_momentum_N;
        draft_drag_N_   = ao.draft_drag_N;
    }
    prof_add(ProfileStage::AeroTruth, 1);

    
// Time-integrated suppression: delivered agent -> exposure -> knockdown -> reduced fire intensity.
    {
        const double delivered_total_kgps_base =
            std::max(0.0, agent_mdot_kgps_ * clamp01(hit_efficiency_0_1_))
            * std::max(0.0, control_inputs_.mdot_scale);
        const auto w = sectorWeightsFromSprayDir(spray_dir_unit_);

        // --------------------
        // Phase 3A: ship-quality geometry realism (stabilized occlusion/shielding + dt-consistent LoA)
        // --------------------
        // Scene anchors mirror the visualizer layout for narrative consistency.
        const Vec3d fire_center_m = hotspot_pos_m_;

        // Sector AABBs (2x2) approximating the Phase 2C fire proxy. Kept constant (no visual coupling).
        const Vec3d fire_half_m{0.35, 0.45, 0.35};
        const Vec3d sector_half_m{fire_half_m.x * 0.48, fire_half_m.y, fire_half_m.z * 0.48};

        const double sx[4] = {-1.0, +1.0, -1.0, +1.0};
        const double sz[4] = {-1.0, -1.0, +1.0, +1.0};

        std::array<Vec3d, kNumSectors_> sector_center_m{};
        std::array<AABBd, kNumSectors_> sector_aabb{};
        for (int i = 0; i < kNumSectors_; ++i) {
            sector_center_m[i] = {fire_center_m.x + sx[i] * sector_half_m.x,
                                  fire_center_m.y,
                                  fire_center_m.z + sz[i] * sector_half_m.z};
            sector_aabb[i] = {sector_center_m[i], sector_half_m};
        }

        // Compute per-sector geometry attenuations.
        // - Occlusion: ray vs scenario-owned obstacle list (up to 2 AABBs), stabilized with hysteresis.
        // - LoA: explicit curved response + dt-consistent smoothing.
        // - Shielding: near sector AABBs shadow far sectors along ray, stabilized with hysteresis.

        // Convert stabilization seconds into deterministic step counts from dt.
        const int enter_steps = std::max(1, (int)std::ceil(hysteresis_enter_s_ / dt));
        const int exit_steps  = std::max(1, (int)std::ceil(hysteresis_exit_s_  / dt));

        sector_net_delivered_mdot_kgps_.fill(0.0);
    sector_utilization_U_0_1_.fill(0.0);
    sector_effective_exposure_kg_.fill(0.0);
    sector_EC50_adj_kg_.fill(0.0);
    utilization_U_0_1_ = 0.0;
    effective_exposure_kg_ = 0.0;
    EC50_adj_kg_ = 0.0;
    telemetry_crc_u32_ = 0;
        sector_raw_delivered_mdot_kgps_.fill(0.0);
        sector_shield_0_1_.fill(1.0);

        for (int j = 0; j < kNumSectors_; ++j) {
            double raw = delivered_total_kgps_base * w[j];
            // Autonomy seam: deterministic sector enable mask (default enables all sectors).
            if (((control_inputs_.sector_enable_mask_u32 >> j) & 1u) == 0u) raw = 0.0;
            sector_raw_delivered_mdot_kgps_[j] = std::max(0.0, raw);

            // --- Occlusion (raw) ---
            bool blocked_raw = false;
            {
                // Phase 3CA.1 scalability: deterministic multi-ray occlusion sampling.
                // Baseline equivalence: rays_per_sector_ defaults to 1 => identical to Phase 3B.1 behavior.
                const int rays = std::clamp(rays_per_sector_, 1, 9);
                static const std::array<Vec3d, 9> kRayOffsetsNorm = {{
                    Vec3d{0.0, 0.0, 0.0},
                    Vec3d{+0.35, 0.0, 0.0}, Vec3d{-0.35, 0.0, 0.0},
                    Vec3d{0.0, 0.0, +0.35}, Vec3d{0.0, 0.0, -0.35},
                    Vec3d{+0.35, 0.0, +0.35}, Vec3d{-0.35, 0.0, +0.35},
                    Vec3d{+0.35, 0.0, -0.35}, Vec3d{-0.35, 0.0, -0.35}
                }};

                bool any_clear = false;
                for (int r = 0; r < rays; ++r) {
                    const Vec3d target = {
                        sector_center_m[j].x + kRayOffsetsNorm[r].x * sector_half_m.x,
                        sector_center_m[j].y,
                        sector_center_m[j].z + kRayOffsetsNorm[r].z * sector_half_m.z
                    };
                    const Vec3d to = {target.x - nozzle_pos_m_.x, target.y - nozzle_pos_m_.y, target.z - nozzle_pos_m_.z};
                    const double dist = std::sqrt(to.x*to.x + to.y*to.y + to.z*to.z);
                    const Vec3d ray_dir = (dist > kEps) ? Vec3d{to.x/dist, to.y/dist, to.z/dist} : Vec3d{0.0,0.0,0.0};

                    bool blocked_this = false;
                    for (int k = 0; k < num_obstacles_; ++k) {
                        AABBd b = obstacles_[k];
                        b.h.x = std::max(0.0, b.h.x + aabb_pad_m_);
                        b.h.y = std::max(0.0, b.h.y + aabb_pad_m_);
                        b.h.z = std::max(0.0, b.h.z + aabb_pad_m_);

                        double tE = 0.0, tX = 0.0;
                        prof_add(ProfileStage::GeometryAttenuation, 1); // ray-AABB test
                        const bool hit = (dist > kEps) && rayAabbIntersect(nozzle_pos_m_, ray_dir, b, tE, tX);
                        if (hit && tE >= 0.0 && tE < (dist - 1e-4)) {
                            blocked_this = true;
                            break;
                        }
                    }

                    if (!blocked_this) {
                        any_clear = true;
                        break;
                    }
                }

                // Multi-ray policy: sector is considered blocked only if all rays are blocked.
                blocked_raw = !any_clear;
            }
            // --- Occlusion (stabilized hysteresis) ---
            {
                const bool want_blocked = blocked_raw;
                const bool is_blocked   = (occ_stable_0_1_[j] < 0.5);

                if (want_blocked) {
                    occ_enter_count_[j]++;
                    occ_exit_count_[j] = 0;
                    if (!is_blocked && occ_enter_count_[j] >= enter_steps) {
                        occ_stable_0_1_[j] = 0.0;
                    }
                } else {
                    occ_exit_count_[j]++;
                    occ_enter_count_[j] = 0;
                    if (is_blocked && occ_exit_count_[j] >= exit_steps) {
                        occ_stable_0_1_[j] = 1.0;
                    }
                }
                sector_occlusion_0_1_[j] = clamp01(occ_stable_0_1_[j]);
            }

            // --- Line-of-Attack (curve + dt-consistent smoothing) ---
            {
                const Vec3d n = safeNorm({nozzle_pos_m_.x - sector_center_m[j].x,
                                         nozzle_pos_m_.y - sector_center_m[j].y,
                                         nozzle_pos_m_.z - sector_center_m[j].z});
                const double loa_target = lineOfAttackCurved_0_1(spray_dir_unit_, n, loa_min_0_1_, loa_power_);
                loa_smooth_0_1_[j] = firstOrderLag(loa_smooth_0_1_[j], loa_target, dt, loa_smooth_tau_s_);
                sector_line_attack_0_1_[j] = clamp01(loa_smooth_0_1_[j]);
            }

            // --- Shielding (raw) ---
            double shield_raw = 1.0;
            {
                const Vec3d toJ = {sector_center_m[j].x - nozzle_pos_m_.x,
                                   sector_center_m[j].y - nozzle_pos_m_.y,
                                   sector_center_m[j].z - nozzle_pos_m_.z};
                const double distJ = std::sqrt(toJ.x*toJ.x + toJ.y*toJ.y + toJ.z*toJ.z);
                const Vec3d dirJ = (distJ > kEps) ? Vec3d{toJ.x/distJ, toJ.y/distJ, toJ.z/distJ} : Vec3d{0.0,0.0,0.0};

                for (int i = 0; i < kNumSectors_; ++i) {
                    if (i == j) continue;
                    const Vec3d toI = {sector_center_m[i].x - nozzle_pos_m_.x,
                                       sector_center_m[i].y - nozzle_pos_m_.y,
                                       sector_center_m[i].z - nozzle_pos_m_.z};
                    const double distI = std::sqrt(toI.x*toI.x + toI.y*toI.y + toI.z*toI.z);
                    if (distI <= kEps || distI >= distJ) continue; // only nearer can shield farther

                    double tE = 0.0, tX = 0.0;
                    prof_add(ProfileStage::GeometryAttenuation, 1); // ray-AABB test (shielding)
                    if (!rayAabbIntersect(nozzle_pos_m_, dirJ, sector_aabb[i], tE, tX)) continue;
                    if (tE >= 0.0 && tE < (distJ - 1e-4)) {
                        constexpr double kShadowLeak = 0.15;
                        shield_raw = std::min(shield_raw, kShadowLeak);
                    }
                }
            }

            // --- Shielding (stabilized hysteresis) ---
            {
                const bool want_shielded = (shield_raw < 0.999);
                const bool is_shielded   = (shield_stable_0_1_[j] < 0.999);

                if (want_shielded) {
                    shield_enter_count_[j]++;
                    shield_exit_count_[j] = 0;
                    if (!is_shielded && shield_enter_count_[j] >= enter_steps) {
                        shield_stable_0_1_[j] = clamp01(shield_raw);
                    }
                } else {
                    shield_exit_count_[j]++;
                    shield_enter_count_[j] = 0;
                    if (is_shielded && shield_exit_count_[j] >= exit_steps) {
                        shield_stable_0_1_[j] = 1.0;
                    }
                }
                sector_shield_0_1_[j] = clamp01(shield_stable_0_1_[j]);
            }

            const double geom = clamp01(sector_occlusion_0_1_[j])
                              * clamp01(sector_line_attack_0_1_[j])
                              * clamp01(sector_shield_0_1_[j]);

            sector_net_delivered_mdot_kgps_[j] = std::max(0.0, raw * geom);
            sector_delivered_mdot_kgps_[j] = sector_net_delivered_mdot_kgps_[j];
        }

        prof_add(ProfileStage::GeometryAttenuation, 1);

        // Optional recovery: decay exposure when delivery stops (after geometry).
        const double sum_raw_kgps = std::max(0.0,
            sector_raw_delivered_mdot_kgps_[0] + sector_raw_delivered_mdot_kgps_[1] +
            sector_raw_delivered_mdot_kgps_[2] + sector_raw_delivered_mdot_kgps_[3]);
        const double net_delivered_total_kgps =
            std::max(0.0, sector_net_delivered_mdot_kgps_[0] + sector_net_delivered_mdot_kgps_[1]
                           + sector_net_delivered_mdot_kgps_[2] + sector_net_delivered_mdot_kgps_[3]);

        const bool delivering = (net_delivered_total_kgps > 1e-6);

        // Safety rails (deterministic hold timers; exported as warnings in observe()).
        // Raw > 0 but net == 0 for >= 0.5 s -> likely fully blocked.
        if (sum_raw_kgps > 1e-6 && net_delivered_total_kgps <= 1e-9) {
            fully_blocked_hold_s_ += dt;
        } else {
            fully_blocked_hold_s_ = 0.0;
        }

        // Sustained glancing hold: net is flowing but LoA remains low.
        double loa_avg = 0.0;
        for (int i = 0; i < kNumSectors_; ++i) loa_avg += clamp01(sector_line_attack_0_1_[i]);
        loa_avg /= (double)kNumSectors_;
        glancing_hold_s_ = (delivering && loa_avg < 0.35) ? (glancing_hold_s_ + dt) : 0.0;
        const double tau_decay = (enable_recovery_ ? std::max(tau_exposure_decay_s_, 0.1) : 1e30);
        const double decay = std::exp(-dt / tau_decay);

        for (int i = 0; i < kNumSectors_; ++i) {
            if (!delivering && enable_recovery_) {
                sector_exposure_kg_[i] *= decay;
            }

            sector_exposure_kg_[i] += sector_delivered_mdot_kgps_[i] * dt;
            if (!std::isfinite(sector_exposure_kg_[i]) || sector_exposure_kg_[i] < 0.0) sector_exposure_kg_[i] = 0.0;

        // Geometry attenuation complete.
        prof_add(ProfileStage::GeometryAttenuation, 1);

        // Phase 3B: chemical effectiveness layer (deterministic, calibratable)
const double U = utilizationU(sector_exposure_kg_[i], agent_profile_.k_util_1_per_kg);
sector_utilization_U_0_1_[i] = U;

const double eff_exposure = sector_exposure_kg_[i] * U * std::max(0.0, agent_profile_.potency);
sector_effective_exposure_kg_[i] = (std::isfinite(eff_exposure) ? std::max(0.0, eff_exposure) : 0.0);

const double EC50_adj = std::max(kEps, agent_profile_.EC50_kg)
    * std::max(kEps, scenario_factors_.temp_factor)
    * std::max(kEps, scenario_factors_.vent_factor)
    * std::max(kEps, scenario_factors_.fuel_factor);
sector_EC50_adj_kg_[i] = EC50_adj;

sector_knockdown_target_0_1_[i] =
    knockdownTargetFromEffectiveExposure(sector_effective_exposure_kg_[i], EC50_adj, agent_profile_.hill, kEps);
            sector_knockdown_0_1_[i] = firstOrderLag(sector_knockdown_0_1_[i], sector_knockdown_target_0_1_[i], dt, tau_knockdown_rise_s_);
        }

        // Aggregate (for Phase 2B semantics)
        exposure_kg_ = 0.0;
        knockdown_0_1_ = 0.0;
        for (int i = 0; i < kNumSectors_; ++i) {
            exposure_kg_ += sector_exposure_kg_[i];
            knockdown_0_1_ += sector_knockdown_0_1_[i];
        }
        knockdown_0_1_ = clamp01(knockdown_0_1_ / (double)kNumSectors_);
// Phase 3B: aggregate chemical layer (for headline + HRR mapping)
utilization_U_0_1_ = utilizationU(exposure_kg_, agent_profile_.k_util_1_per_kg);
effective_exposure_kg_ = (std::isfinite(exposure_kg_) ? std::max(0.0, exposure_kg_) : 0.0)
    * utilization_U_0_1_ * std::max(0.0, agent_profile_.potency);

EC50_adj_kg_ = std::max(kEps, agent_profile_.EC50_kg)
    * std::max(kEps, scenario_factors_.temp_factor)
    * std::max(kEps, scenario_factors_.vent_factor)
    * std::max(kEps, scenario_factors_.fuel_factor);

knockdown_target_0_1_ =
    knockdownTargetFromEffectiveExposure(effective_exposure_kg_, EC50_adj_kg_, agent_profile_.hill, kEps);

        // Regime classification based on net delivered mass flow (after geometry).
        if (net_delivered_total_kgps <= 1e-6) suppression_regime_ = SuppressionRegime::None;
        else if (net_delivered_total_kgps < 0.02) suppression_regime_ = SuppressionRegime::Ineffective;
        else if (net_delivered_total_kgps < 0.06) suppression_regime_ = SuppressionRegime::Marginal;
        else if (net_delivered_total_kgps < 0.12) suppression_regime_ = SuppressionRegime::Effective;
        else suppression_regime_ = SuppressionRegime::Overkill;
    }
    prof_add(ProfileStage::ExposureKnockdown, 1);

    // Geometry -> exposure -> knockdown -> HRR mapping is complete for this step.
    prof_add(ProfileStage::ExposureKnockdown, 1);

    // --------------------
    // Li-ion runaway
    // --------------------
    liion_.step(dt, reactor_.temperatureK(), liionHeat_W_, liionVent_kgps_);

    if (liionVent_kgps_ > 0.0) {
        const double m = liionVent_kgps_ * dt;

        const auto& sp = reactor_.species();
        const double Minert = sp[idx_.iINERT].molarMass_kg_per_mol;

        if (Minert > kEps) {
            reactor_.addMoles(idx_.iINERT, m / Minert);
        }
    }

    // External cooling "felt" by the reactor is reduced by an internal heat source.
    const double effectiveExternalCooling_W = cooling_W - liionHeat_W_;

    // --------------------
    // Reactor step
    // --------------------
    double combustionHRR_raw_W = 0.0;

// Phase 2B/2C causal coupling: knockdown affects the *thermal* coupling deterministically.
const double combustionHeatMult_0_1 = 1.0 - clamp01(knockdown_0_1_);

// Reactor returns RAW combustion HRR (pre-multiplier). The multiplier is applied to the thermal state internally.
// Post-ignition, provide a small kinetics assist so combustion can start deterministically at ambient,
// without spoofing telemetry: HRR remains chemistry-derived and is applied consistently to the thermal state.
const double ignitionTempFloor_K = ignited_ ? 600.0 : 0.0;
reactor_.step(dt, inhib_kgm3_, effectiveExternalCooling_W, combustionHeatMult_0_1, combustionHRR_raw_W, ignitionTempFloor_K);
prof_add(ProfileStage::ReactorStep, 1);

// Effective combustion HRR is raw scaled by multiplier (telemetry aligned with thermal coupling).
double combustionHRR_W = combustionHRR_raw_W * combustionHeatMult_0_1;
if (!std::isfinite(combustionHRR_W) || combustionHRR_W < 0.0) combustionHRR_W = 0.0;

    raw_HRR_W_ = std::max(0.0, combustionHRR_raw_W + liionHeat_W_);
    effective_HRR_W_ = std::max(0.0, combustionHRR_W + liionHeat_W_);


    // Net HRR = combustion + li-ion contribution
    lastHRR_W_ = effective_HRR_W_;

    prof_add(ProfileStage::ReactorStep, 1);


// --------------------
// Phase 3B.1 harness sampler (deterministic, sim-time based)
// --------------------
// Note: this must not affect the simulation state.
// Sampling happens at fixed sim-time cadence (telemetry_dt_s_).


    // --------------------
    // Ventilation step
    // --------------------
    vent_.apply(dt, reactor_);
    prof_add(ProfileStage::Ventilation, 1);
    prof_add(ProfileStage::Ventilation, 1);

    // --------------------
    // Recompute end-of-step inert concentration for visualization correctness
    // --------------------
    {
        const double V = reactor_.config().volume_m3;
        inert_kgm3_ = 0.0;

        if (V > kEps) {
            const auto& sp = reactor_.species();
            const auto& n  = reactor_.moles();

            if (idx_.iINERT >= 0
                && idx_.iINERT < static_cast<int>(n.size())
                && idx_.iINERT < static_cast<int>(sp.size())) {

                const double M = sp[idx_.iINERT].molarMass_kg_per_mol;
                const double nInert = n[idx_.iINERT];

                if (M > kEps && std::isfinite(nInert) && nInert > 0.0) {
                    inert_kgm3_ = (nInert * M) / V;
                }
            }
        }
    }

    // --------------------
    // Phase 3B.1: deterministic telemetry sampling (schema v1)
    // --------------------
    {
        auto fnv_add_f32 = [&](std::uint32_t h, float v) {
            std::uint32_t bits = 0;
            std::memcpy(&bits, &v, sizeof(v));
            return fnv1a32_update(h, &bits, sizeof(bits));
        };

        auto push_sample = [&](double sample_t_s) {
            TelemetrySampleV1 s{};
            s.t_s = static_cast<float>(sample_t_s);

            const double sum_raw =
                sector_raw_delivered_mdot_kgps_[0] + sector_raw_delivered_mdot_kgps_[1] +
                sector_raw_delivered_mdot_kgps_[2] + sector_raw_delivered_mdot_kgps_[3];
            const double sum_net =
                sector_net_delivered_mdot_kgps_[0] + sector_net_delivered_mdot_kgps_[1] +
                sector_net_delivered_mdot_kgps_[2] + sector_net_delivered_mdot_kgps_[3];

            s.raw_mdot_kgps = static_cast<float>(std::max(0.0, sum_raw));
            s.net_mdot_kgps = static_cast<float>(std::max(0.0, sum_net));
            s.exposure_kg = static_cast<float>(std::max(0.0, exposure_kg_));
            s.effective_exposure_kg = static_cast<float>(std::max(0.0, effective_exposure_kg_));
            s.KD_target_0_1 = static_cast<float>(clamp01(knockdown_target_0_1_));
            s.KD_actual_0_1 = static_cast<float>(clamp01(knockdown_0_1_));
            s.HRR_kW = static_cast<float>(std::max(0.0, effective_HRR_W_) * 0.001);

            telemetry_rb_[telemetry_head_] = s;
            telemetry_head_ = (telemetry_head_ + 1) % kTelemetryCapacity_;
            if (telemetry_count_ < kTelemetryCapacity_) telemetry_count_++;

            std::uint32_t events = 0u;
            const bool occluded_any =
                (sector_occlusion_0_1_[0] < 0.5) || (sector_occlusion_0_1_[1] < 0.5) ||
                (sector_occlusion_0_1_[2] < 0.5) || (sector_occlusion_0_1_[3] < 0.5);
            if (!prev_occluded_any_ && occluded_any) events |= Event_OcclusionEnter;
            if (prev_occluded_any_ && !occluded_any) events |= Event_OcclusionExit;
            prev_occluded_any_ = occluded_any;

            const double kd_now = clamp01(knockdown_0_1_);
            const bool kd_ge_0p5 = (kd_now >= 0.5);
            const bool kd_ge_0p9 = (kd_now >= 0.9);
            if (!prev_kd_ge_0p5_ && kd_ge_0p5) events |= Event_KD_ge_0p5;
            if (!prev_kd_ge_0p9_ && kd_ge_0p9) events |= Event_KD_ge_0p9;
            prev_kd_ge_0p5_ = kd_ge_0p5;
            prev_kd_ge_0p9_ = kd_ge_0p9;

            const bool hrr_below_100kW = (std::max(0.0, effective_HRR_W_) < 100000.0);
            if (!prev_hrr_below_100kW_ && hrr_below_100kW) events |= Event_HRR_below_100kW;
            prev_hrr_below_100kW_ = hrr_below_100kW;

            if (sum_net > sum_raw + 1e-9) events |= Warn_NetMdot_gt_RawMdot;
            if (std::isfinite(prev_exposure_kg_) && exposure_kg_ + 1e-12 < prev_exposure_kg_) events |= Warn_ExposureDecreased;
            if (std::isfinite(prev_effective_exposure_kg_) && std::isfinite(prev_kd_target_0_1_)) {
                if (effective_exposure_kg_ > prev_effective_exposure_kg_ + 1e-9
                    && knockdown_target_0_1_ + 1e-9 < prev_kd_target_0_1_) {
                    events |= Warn_KD_NonMonotonic;
                }
            }
            prev_exposure_kg_ = exposure_kg_;
            prev_effective_exposure_kg_ = effective_exposure_kg_;
            prev_kd_target_0_1_ = knockdown_target_0_1_;
            latest_events_bits_ |= events;

            // Phase 3CA: telemetry schema v2 (derived monitoring, never hashed)
            if (telemetry_v2_enabled_) {
                TelemetrySampleV2 s2{};
                s2.t_s = s.t_s;
                s2.raw_mdot_kgps = s.raw_mdot_kgps;
                s2.net_mdot_kgps = s.net_mdot_kgps;
                s2.exposure_kg = s.exposure_kg;
                s2.effective_exposure_kg = s.effective_exposure_kg;
                s2.KD_target_0_1 = s.KD_target_0_1;
                s2.KD_actual_0_1 = s.KD_actual_0_1;
                s2.HRR_kW = s.HRR_kW;
                s2.events_u32 = events;

                double occ_sum = 0.0, occ_min = 1.0;
                double loa_sum = 0.0, loa_min = 1.0;
                double sh_sum = 0.0, sh_min = 1.0;
                std::uint32_t blocked = 0u;
                for (int i = 0; i < 4; ++i) {
                    const double occ = clamp01(sector_occlusion_0_1_[i]);
                    const double loa = clamp01(sector_line_attack_0_1_[i]);
                    const double sh  = clamp01(sector_shield_0_1_[i]);
                    occ_sum += occ;
                    loa_sum += loa;
                    sh_sum  += sh;
                    occ_min = std::min(occ_min, occ);
                    loa_min = std::min(loa_min, loa);
                    sh_min  = std::min(sh_min, sh);
                    if (occ < 0.5) blocked++;
                }
                s2.occ_avg_0_1 = static_cast<float>(occ_sum * 0.25);
                s2.occ_min_0_1 = static_cast<float>(occ_min);
                s2.loa_avg_0_1 = static_cast<float>(loa_sum * 0.25);
                s2.loa_min_0_1 = static_cast<float>(loa_min);
                s2.shield_avg_0_1 = static_cast<float>(sh_sum * 0.25);
                s2.shield_min_0_1 = static_cast<float>(sh_min);
                s2.blocked_sector_count_u32 = blocked;

                telemetry_v2_rb_[telemetry_v2_head_] = s2;
                telemetry_v2_head_ = (telemetry_v2_head_ + 1) % kTelemetryV2Capacity_;
                if (telemetry_v2_count_ < kTelemetryV2Capacity_) telemetry_v2_count_++;
            }

            if (verification_mode_ || calibration_mode_) {
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.t_s);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.raw_mdot_kgps);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.net_mdot_kgps);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.exposure_kg);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.effective_exposure_kg);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.KD_target_0_1);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.KD_actual_0_1);
                run_signatures_.telemetry_crc_u32 = crc32_add_f32(run_signatures_.telemetry_crc_u32, s.HRR_kW);

                std::uint32_t h = run_signatures_.state_digest_u32;
                if (h == 0) h = fnv1a32_begin();
                h = fnv_add_f32(h, s.t_s);
                h = fnv_add_f32(h, s.exposure_kg);
                h = fnv_add_f32(h, s.effective_exposure_kg);
                h = fnv_add_f32(h, s.KD_target_0_1);
                h = fnv_add_f32(h, s.KD_actual_0_1);
                h = fnv_add_f32(h, s.HRR_kW);
                for (int i = 0; i < 4; ++i) {
                    h = fnv_add_f32(h, static_cast<float>(clamp01(sector_occlusion_0_1_[i])));
                    h = fnv_add_f32(h, static_cast<float>(clamp01(sector_line_attack_0_1_[i])));
                    h = fnv_add_f32(h, static_cast<float>(clamp01(sector_shield_0_1_[i])));
                    h = fnv_add_f32(h, static_cast<float>(std::max(0.0, sector_exposure_kg_[i])));
                    h = fnv_add_f32(h, static_cast<float>(std::max(0.0, sector_effective_exposure_kg_[i])));
                    h = fnv_add_f32(h, static_cast<float>(clamp01(sector_knockdown_target_0_1_[i])));
                    h = fnv_add_f32(h, static_cast<float>(clamp01(sector_knockdown_0_1_[i])));
                }
                run_signatures_.state_digest_u32 = h;
            }
        };

        if (!std::isfinite(telemetry_dt_s_) || telemetry_dt_s_ < 1e-6) telemetry_dt_s_ = 0.10;
        while (scenario_time_s_ + 1e-12 >= telemetry_next_t_s_) {
            push_sample(telemetry_next_t_s_);
            telemetry_next_t_s_ += telemetry_dt_s_;
        }
    }

    prof_add(ProfileStage::TelemetrySample, 1);
    if (profiling_enabled_) {
        last_profile_.t_s = static_cast<float>(scenario_time_s_);
    }


    // --------------------
    // Termination logic
    // --------------------
    const double T_C = reactor_.temperatureK() - 273.15;

    // Original strict conclude condition (fuel nearly exhausted)
    if (!verification_mode_ && ignited_
        && lastHRR_W_ < kConclude_HRR_W
        && T_C < kConclude_T_C
        && fuelSolid_kg_ < kConclude_fuel_kg) {
        concluded_ = true;
        return;
    }

    // UI-friendly safe-hold termination:
    // require "safe" conditions to persist continuously for kConcludeSafeHold_s.
    const bool suppressionLikelyActive =
        (agent_mdot_kgps_ > 0.0) || (inhib_kgm3_ > 0.0) || (inert_kgm3_ > 0.0);

    const bool safeNow =
        ignited_
        && suppressionLikelyActive
        && (lastHRR_W_ < kConcludeSafe_HRR_W)
        && (T_C < kConcludeSafe_T_C);

    if (safeNow) {
        safeHold_s_ += dt;
    } else {
        safeHold_s_ = 0.0;
    }

    if (!verification_mode_ && safeHold_s_ >= kConcludeSafeHold_s) {
        concluded_ = true;
        return;
    }
}

Observation Simulation::observe() const {
    Observation o;
    o.T_K   = reactor_.temperatureK();
    o.HRR_W = lastHRR_W_;

    // Phase 2B/2C: raw/effective HRR bookkeeping + suppression telemetry
    o.raw_HRR_W = std::isfinite(raw_HRR_W_) ? raw_HRR_W_ : 0.0;
    o.effective_HRR_W = std::isfinite(effective_HRR_W_) ? effective_HRR_W_ : o.HRR_W;

    const double delivered_total_kgps = std::max(0.0, agent_mdot_kgps_ * clamp01(hit_efficiency_0_1_));
    o.delivered_mdot_kgps = std::isfinite(delivered_total_kgps) ? delivered_total_kgps : 0.0;
    // Phase 3: net delivered is the post-geometry sum of sector net deliveries.
    const double net_delivered_total_kgps = std::max(0.0,
        sector_net_delivered_mdot_kgps_[0] + sector_net_delivered_mdot_kgps_[1] +
        sector_net_delivered_mdot_kgps_[2] + sector_net_delivered_mdot_kgps_[3]);
    o.net_delivered_mdot_kgps = std::isfinite(net_delivered_total_kgps) ? net_delivered_total_kgps : 0.0;
o.exposure_kg = (std::isfinite(exposure_kg_) ? std::max(0.0, exposure_kg_) : 0.0);

// Phase 3B: chemical layer aggregates
o.utilization_U_0_1 = clamp01(utilization_U_0_1_);
o.effective_exposure_kg = (std::isfinite(effective_exposure_kg_) ? std::max(0.0, effective_exposure_kg_) : 0.0);
o.EC50_adj_kg = (std::isfinite(EC50_adj_kg_) ? std::max(0.0, EC50_adj_kg_) : 0.0);

o.agent_type = static_cast<int>(agent_type_);
o.agent_EC50_kg = agent_profile_.EC50_kg;
o.agent_hill = agent_profile_.hill;
o.agent_k_util_1_per_kg = agent_profile_.k_util_1_per_kg;
o.agent_potency = agent_profile_.potency;

o.temp_factor = scenario_factors_.temp_factor;
o.vent_factor = scenario_factors_.vent_factor;
o.fuel_factor = scenario_factors_.fuel_factor;

o.run_param_hash_u32 = run_signatures_.run_param_hash_u32;
o.telemetry_crc_u32  = run_signatures_.telemetry_crc_u32;
o.calibration_mode   = calibration_mode_;

    o.knockdown_0_1 = clamp01(knockdown_0_1_);
    o.suppression_regime = static_cast<int>(suppression_regime_);

    // Phase 3A aggregate telemetry (computed from stabilized per-sector state)
    double occ_sum = 0.0;
    double occ_max = 0.0;
    double loa_sum = 0.0;
    double loa_min = 1.0;
    double sum_raw = 0.0;
    double sum_net = 0.0;
    int blocked = 0, shielded = 0, glancing = 0, direct = 0;

    for (int i = 0; i < Observation::kNumSuppressionSectors; ++i) {
        o.sector_delivered_mdot_kgps[i] = (std::isfinite(sector_delivered_mdot_kgps_[i]) ? std::max(0.0, sector_delivered_mdot_kgps_[i]) : 0.0);
        o.sector_exposure_kg[i] = (std::isfinite(sector_exposure_kg_[i]) ? std::max(0.0, sector_exposure_kg_[i]) : 0.0);
        o.sector_knockdown_0_1[i] = clamp01(sector_knockdown_0_1_[i]);
        const double occ = clamp01(sector_occlusion_0_1_[i]);
        const double loa = clamp01(sector_line_attack_0_1_[i]);
        const double net = (std::isfinite(sector_net_delivered_mdot_kgps_[i]) ? std::max(0.0, sector_net_delivered_mdot_kgps_[i]) : 0.0);
        const double raw = (std::isfinite(sector_raw_delivered_mdot_kgps_[i]) ? std::max(0.0, sector_raw_delivered_mdot_kgps_[i]) : 0.0);
        const double sh  = clamp01(sector_shield_0_1_[i]);

        o.sector_occlusion_0_1[i] = occ;
        o.sector_line_attack_0_1[i] = loa;
        o.sector_net_delivered_mdot_kgps[i] = net;

        // Phase 3A auditability fields
        o.sector_shield_0_1[i] = sh;
        o.sector_raw_delivered_mdot_kgps[i] = raw;
        o.sector_knockdown_target_0_1[i] = clamp01(sector_knockdown_target_0_1_[i]);
// Phase 3B: chemical telemetry
o.sector_utilization_U_0_1[i] = clamp01(sector_utilization_U_0_1_[i]);
o.sector_effective_exposure_kg[i] = (std::isfinite(sector_effective_exposure_kg_[i]) ? std::max(0.0, sector_effective_exposure_kg_[i]) : 0.0);
o.sector_EC50_adj_kg[i] = (std::isfinite(sector_EC50_adj_kg_[i]) ? std::max(0.0, sector_EC50_adj_kg_[i]) : 0.0);


        occ_sum += occ;
        occ_max = std::max(occ_max, occ);
        loa_sum += loa;
        loa_min = std::min(loa_min, loa);
        sum_raw += raw;
        sum_net += net;

        // Only classify sectors that are receiving non-trivial raw delivery (story-first).
        if (raw > 1e-6) {
            if (occ < 0.5) blocked++;
            else if (sh < 0.999) shielded++;
            else if (loa < 0.75) glancing++;
            else direct++;
        }
    }

    o.occ_avg_0_1 = clamp01(occ_sum / (double)Observation::kNumSuppressionSectors);
    o.occ_max_0_1 = clamp01(occ_max);
    o.loa_avg_0_1 = clamp01(loa_sum / (double)Observation::kNumSuppressionSectors);
    o.loa_min_0_1 = clamp01(loa_min);
    o.sum_raw_mdot_kgps = std::max(0.0, sum_raw);
    o.sum_net_mdot_kgps = std::max(0.0, sum_net);
    o.blocked_sector_count = blocked;
    o.shielded_sector_count = shielded;
    o.glancing_sector_count = glancing;
    o.direct_sector_count = direct;

    // Headline classification (priority order: Blocked > Shielded > Glancing > Direct)
    if (blocked > 0) o.headline_state = 2;
    else if (shielded > 0) o.headline_state = 3;
    else if (glancing > 0) o.headline_state = 1;
    else o.headline_state = 0;

    o.warn_fully_blocked = (fully_blocked_hold_s_ >= 0.50);
    o.warn_glancing_hold = (glancing_hold_s_ >= 2.00);


    o.O2_volpct  = reactor_.gasMoleFraction(idx_.iO2)  * 100.0;
    o.CO2_volpct = reactor_.gasMoleFraction(idx_.iCO2) * 100.0;
    o.H2O_volpct = reactor_.gasMoleFraction(idx_.iH2O) * 100.0;

    o.fuel_kg        = fuelSolid_kg_;
    o.inhibitor_kgm3 = inhib_kgm3_;
    o.inert_kgm3     = inert_kgm3_;

    o.ACH             = vent_.config().ACH;
    o.agent_mdot_kgps = agent_mdot_kgps_;

    // Phase 2A truth telemetry
    o.vfep_rpm            = std::isfinite(vfep_rpm_) ? vfep_rpm_ : 0.0;
    o.hit_efficiency_0_1  = std::isfinite(hit_efficiency_0_1_) ? std::clamp(hit_efficiency_0_1_, 0.0, 1.0) : 0.0;

    o.spray_dir_unit_x = (std::isfinite(spray_dir_unit_.x) ? spray_dir_unit_.x : 0.0);
    o.spray_dir_unit_y = (std::isfinite(spray_dir_unit_.y) ? spray_dir_unit_.y : 0.0);
    o.spray_dir_unit_z = (std::isfinite(spray_dir_unit_.z) ? spray_dir_unit_.z : 0.0);

    o.draft_vel_mps_x = (std::isfinite(draft_vel_mps_.x) ? draft_vel_mps_.x : 0.0);
    o.draft_vel_mps_y = (std::isfinite(draft_vel_mps_.y) ? draft_vel_mps_.y : 0.0);
    o.draft_vel_mps_z = (std::isfinite(draft_vel_mps_.z) ? draft_vel_mps_.z : 0.0);

    o.jet_momentum_N = std::isfinite(jet_momentum_N_) ? jet_momentum_N_ : 0.0;
    o.draft_drag_N   = std::isfinite(draft_drag_N_)   ? draft_drag_N_   : 0.0;

    // Reward: temperature safety minus agent usage penalty.
    const double T_C = o.T_K - 273.15;
    const double safety = (T_C < kRewardSafe_T_C) ? kRewardSafeBonus : kRewardUnsafePenalty;
    const double waste  = kRewardWasteCoeff * o.agent_mdot_kgps;

    o.reward = safety - waste;

    // Fire / hotspot truth (meters, world frame)
    o.hotspot_pos_m_x = hotspot_pos_m_.x;
    o.hotspot_pos_m_y = hotspot_pos_m_.y;
    o.hotspot_pos_m_z = hotspot_pos_m_.z;

    return o;
}

} // namespace vfep