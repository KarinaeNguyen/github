#pragma once

#include <vector>
#include <array>
#include <cstdint>

#include "Chemistry.h"
#include "Reactor.h"
#include "Ventilation.h"
#include "Suppression.h"
#include "LiIonRunaway.h"
#include "Aerodynamics.h"

namespace vfep {

// ============================================================
// Phase 3B.1: Deterministic safety harness (verification + schema contracts)
//
// Design intent:
// - Add verification infrastructure and stable interfaces.
// - Do NOT change Phase 3B chemical layer math or causal chain.
// - No new dependencies; no architecture rewrites.
// ============================================================

enum class VerificationTestId : std::uint32_t {
    V0 = 0, // Chemistry-dominated (no obstacles, direct LoA)
    V1 = 1, // Geometry-dominated (blocked LoA for at least 1–2 sectors)
    V2 = 2, // Hysteresis edge (threshold crossings via deterministic nozzle sweep)
};

struct TelemetrySampleV1 {
    // Fixed order, explicit units.
    float t_s = 0.0f;
    float raw_mdot_kgps = 0.0f;
    float net_mdot_kgps = 0.0f;
    float exposure_kg = 0.0f;
    float effective_exposure_kg = 0.0f;
    float KD_target_0_1 = 0.0f;
    float KD_actual_0_1 = 0.0f;
    float HRR_kW = 0.0f; // effective HRR in kW
};

// ============================================================
// Phase 3CA: Telemetry schema v2 (additive; v1 is immutable)
//
// Rules:
// - Telemetry v1 sampling, CRC, and ring buffer remain unchanged.
// - v2 is strictly additional monitoring (derived, deterministic).
// - v2 must not affect simulation state, signatures, or verification results.
// ============================================================
struct TelemetrySampleV2 {
    // Matches v1 "core" fields for easy alignment + adds derived monitoring.
    float t_s = 0.0f;
    float raw_mdot_kgps = 0.0f;
    float net_mdot_kgps = 0.0f;
    float exposure_kg = 0.0f;
    float effective_exposure_kg = 0.0f;
    float KD_target_0_1 = 0.0f;
    float KD_actual_0_1 = 0.0f;
    float HRR_kW = 0.0f;

    // Deterministic derived event/warning bits (same semantics as TelemetryEventBits).
    std::uint32_t events_u32 = 0u;

    // Deterministic aggregates (audit-friendly, no new physics).
    float occ_avg_0_1 = 1.0f;
    float occ_min_0_1 = 1.0f;
    float loa_avg_0_1 = 0.0f;
    float loa_min_0_1 = 0.0f;
    float shield_avg_0_1 = 1.0f;
    float shield_min_0_1 = 1.0f;
    std::uint32_t blocked_sector_count_u32 = 0u;
};

// ============================================================
// Phase 3CA: Autonomy-ready seams (no embedded AI; deterministic I/O contracts)
// ============================================================
struct ControlInputsV1 {
    std::uint32_t version_u32 = 1;
    std::uint32_t size_bytes_u32 = sizeof(ControlInputsV1);

    // Bitmask enabling/disabling delivery (1 = enabled). Bit i maps to sector i.
    // For Phase 3CA, the default mask enables all legacy sectors.
    std::uint32_t sector_enable_mask_u32 = 0xFu;

    // Deterministic scalar applied to commanded mdot (dimensionless).
    // 1.0 = passthrough.
    double mdot_scale = 1.0;
};

// ============================================================
// Phase 3CA: Deterministic profiling hooks (zero impact on results)
// Profiling output is explicitly excluded from verification hashes.
// ============================================================
enum class ProfileStage : std::uint32_t {
    ScenarioSweep = 0,
    PyrolysisFuel,
    SuppressionApply,
    AeroTruth,
    GeometryAttenuation,
    ExposureKnockdown,
    ReactorStep,
    Ventilation,
    TelemetrySample,
    Count
};

struct ProfileSampleV1 {
    // Phase 3CA.1: Deterministic profiling counters (work units), not wall-clock time.
    // These counters are observational-only and explicitly excluded from verification hashes.
    std::uint32_t version_u32 = 2;
    std::uint32_t size_bytes_u32 = sizeof(ProfileSampleV1);
    float t_s = 0.0f;
    // Deterministic work units per stage (e.g., ray-AABB tests, loop iterations).
    std::array<std::uint64_t, static_cast<std::size_t>(ProfileStage::Count)> stage_work{{0}};
};

struct RunSignatures {
    std::uint32_t run_param_hash_u32 = 0; // FNV-1a32 over effective parameters
    std::uint32_t telemetry_crc_u32  = 0; // CRC32 over sampled telemetry stream (TelemetrySampleV1)
    std::uint32_t state_digest_u32   = 0; // FNV-1a32 over compact state snapshots
};

// Event / warning bitmask (developer-visible; must not affect simulation results).
enum TelemetryEventBits : std::uint32_t {
    Event_None               = 0u,
    Event_OcclusionEnter     = 1u << 0,
    Event_OcclusionExit      = 1u << 1,
    Event_KD_ge_0p5          = 1u << 2,
    Event_KD_ge_0p9          = 1u << 3,
    Event_HRR_below_100kW    = 1u << 4,
    // Invariant checks (warnings)
    Warn_NetMdot_gt_RawMdot  = 1u << 16,
    Warn_ExposureDecreased   = 1u << 17,
    Warn_KD_NonMonotonic     = 1u << 18,
};

// Versioned, hashable, auditable contracts.
// NOTE: These blocks are intentionally minimal and only cover parameters owned by Simulation.
struct SimConfigV1 {
    std::uint32_t version_u32 = 1;
    std::uint32_t size_bytes_u32 = sizeof(SimConfigV1);
    std::uint32_t fnv_hash_u32 = 0;

    // Deterministic stepping + telemetry cadence
    double fixed_dt_s = 0.0;
    double telemetry_dt_s = 0.0;

    // Tunables (explicit units)
    double exposure_half_kg = 0.0;
    double exposure_hill_n = 0.0;
    double tau_knockdown_rise_s = 0.0;
    double tau_knockdown_fall_s = 0.0;
    std::uint32_t enable_recovery_u32 = 0;
    double tau_exposure_decay_s = 0.0;
};

struct AgentConfigV1 {
    std::uint32_t version_u32 = 1;
    std::uint32_t size_bytes_u32 = sizeof(AgentConfigV1);
    std::uint32_t fnv_hash_u32 = 0;

    std::int32_t agent_type_i32 = 0; // AgentType
    double EC50_kg = 0.0;
    double hill = 0.0;
    double k_util_1_per_kg = 0.0;
    double potency = 1.0;
};

struct ScenarioConfigV1 {
    std::uint32_t version_u32 = 1;
    std::uint32_t size_bytes_u32 = sizeof(ScenarioConfigV1);
    std::uint32_t fnv_hash_u32 = 0;

    std::int32_t demo_scenario_i32 = 0; // DemoScenario

    // Nozzle pose (meters, unit vector)
    Vec3d nozzle_pos_m{0.0, 0.0, 0.0};
    Vec3d nozzle_dir_unit{0.0, 0.0, 1.0};

    // Deterministic sweep (if enabled)
    std::uint32_t nozzle_sweep_enabled_u32 = 0;
    double sweep_amp_deg = 0.0;
    double sweep_freq_hz = 0.0;

    // Discharge + ventilation
    double mdot_cmd_kgps = 0.0;
    double ACH_1_per_h = 0.0;

    // Scenario static factors applied to EC50 (dimensionless)
    double temp_factor = 1.0;
    double vent_factor = 1.0;
    double fuel_factor = 1.0;
};

// ---- Phase 3A deterministic demo scenarios ----
enum class DemoScenario : int {
    DirectVsGlance = 0,
    OcclusionWall  = 1,
    ShieldingStack = 2,
    Mixed          = 3,
};

// Minimal axis-aligned box proxy (double precision, simulation-owned).
struct AABBd { Vec3d c; Vec3d h; };


// ---- Phase 3B chemically calibrated effectiveness layer ----
enum class AgentType : int {
    CleanAgent   = 0,
    DryChemical  = 1,
    CO2          = 2,
};

struct AgentProfile {
    // Exposure at 50% knockdown (kg)
    double EC50_kg = 1.0;
    // Hill slope (dimensionless)
    double hill = 1.0;
    // Utilization ramp rate (1/kg)
    double k_util_1_per_kg = 1.0;
    // Optional potency scalar (dimensionless). 1.0 = baseline.
    double potency = 1.0;
};

// Deterministic, static scenario factors (pure scalars). These adjust EC50, not exposure.
struct ScenarioFactors {
    double temp_factor = 1.0; // dimensionless
    double vent_factor = 1.0; // dimensionless
    double fuel_factor = 1.0; // dimensionless
};


// ---- Phase 2B/2C suppression regime ----
enum class SuppressionRegime : int {
    None = 0,
    Ineffective = 1,
    Marginal = 2,
    Effective = 3,
    Overkill = 4,
};

struct Observation {
    double T_K = 295.15;
    double HRR_W = 0.0;
    double O2_volpct = 20.95;
    double CO2_volpct = 0.042;
    double H2O_volpct = 1.0;
    double fuel_kg = 0.0;

    double inhibitor_kgm3 = 0.0;
    double inert_kgm3 = 0.0;

    double ACH = 0.0;
    double agent_mdot_kgps = 0.0;

    // Fire / hotspot truth (meters, world frame)
    double hotspot_pos_m_x = 0.0;
    double hotspot_pos_m_y = 0.0;
    double hotspot_pos_m_z = 0.0;

    // ---- Phase 2A truth telemetry (simulation-produced) ----
    double vfep_rpm = 0.0;                 // actuator truth
    double hit_efficiency_0_1 = 0.0;       // computed in sim/aero module
    double spray_dir_unit_x = 0.0;
    double spray_dir_unit_y = 0.0;
    double spray_dir_unit_z = 0.0;
    double draft_vel_mps_x = 0.0;
    double draft_vel_mps_y = 0.0;
    double draft_vel_mps_z = 0.0;

    // Optional debug telemetry
    double jet_momentum_N = 0.0;
    double draft_drag_N   = 0.0;

    // ---- Phase 2B time-integrated suppression telemetry ----
    double delivered_mdot_kgps = 0.0;
    // Phase 3: net delivered after geometry effects (sum of sector_net_delivered_mdot_kgps)
    double net_delivered_mdot_kgps = 0.0;
    double exposure_kg = 0.0;
    double knockdown_0_1 = 0.0;
    double raw_HRR_W = 0.0;
    double effective_HRR_W = 0.0;
    int suppression_regime = 0;

    // ---- Phase 2C spatial suppression (4 sectors) ----
    static constexpr int kNumSuppressionSectors = 4;
    std::array<double, kNumSuppressionSectors> sector_delivered_mdot_kgps{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSuppressionSectors> sector_exposure_kg{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSuppressionSectors> sector_knockdown_0_1{{0.0,0.0,0.0,0.0}};

    // ---- Phase 3 geometry suppression telemetry (per sector) ----
    // 1.0 = unobstructed; 0.0 = fully occluded by geometry
    std::array<double, kNumSuppressionSectors> sector_occlusion_0_1{{1.0,1.0,1.0,1.0}};
    // 0..1 scalar capturing line-of-attack effectiveness (direct vs glancing)
    std::array<double, kNumSuppressionSectors> sector_line_attack_0_1{{0.0,0.0,0.0,0.0}};
    // Net delivered mass flow after geometry effects (occlusion, line-of-attack, shielding)
    std::array<double, kNumSuppressionSectors> sector_net_delivered_mdot_kgps{{0.0,0.0,0.0,0.0}};

    // ---- Phase 3A ship-quality audit telemetry (per sector) ----
    // Shield scalar applied causally to net delivered mdot.
    std::array<double, kNumSuppressionSectors> sector_shield_0_1{{1.0,1.0,1.0,1.0}};
    // Raw delivered per sector prior to geometry (occlusion/LoA/shield).
    std::array<double, kNumSuppressionSectors> sector_raw_delivered_mdot_kgps{{0.0,0.0,0.0,0.0}};
    // Per-sector knockdown target (Hill), and actual smoothed state.
    std::array<double, kNumSuppressionSectors> sector_knockdown_target_0_1{{0.0,0.0,0.0,0.0}};


// ---- Phase 3B chemical effectiveness telemetry (per sector) ----
// Utilization U(x) = 1 - exp(-k_util * exposure_kg)
std::array<double, kNumSuppressionSectors> sector_utilization_U_0_1{{0.0,0.0,0.0,0.0}};
// Effective exposure = exposure_kg * U(exposure_kg)
std::array<double, kNumSuppressionSectors> sector_effective_exposure_kg{{0.0,0.0,0.0,0.0}};
// Scenario-adjusted EC50 used for the dose-response curve (kg)
std::array<double, kNumSuppressionSectors> sector_EC50_adj_kg{{0.0,0.0,0.0,0.0}};

// Aggregate chemical telemetry
double utilization_U_0_1 = 0.0;
double effective_exposure_kg = 0.0;
double EC50_adj_kg = 0.0;

// Active agent profile (science-ready, explicit units)
int agent_type = 0; // AgentType
double agent_EC50_kg = 0.0;
double agent_hill = 0.0;
double agent_k_util_1_per_kg = 0.0;
double agent_potency = 1.0;

// Scenario static factors applied to EC50 (dimensionless)
double temp_factor = 1.0;
double vent_factor = 1.0;
double fuel_factor = 1.0;

// Calibration signatures (science-ready repeatability)
std::uint32_t run_param_hash_u32 = 0;
std::uint32_t telemetry_crc_u32  = 0;
bool calibration_mode = false;


    // ---- Phase 3A aggregate telemetry (auditability + HUD story) ----
    double occ_avg_0_1 = 1.0;
    double occ_max_0_1 = 1.0;
    double loa_avg_0_1 = 0.0;
    double loa_min_0_1 = 0.0;
    double sum_raw_mdot_kgps = 0.0;
    double sum_net_mdot_kgps = 0.0;
    int blocked_sector_count = 0;
    int shielded_sector_count = 0;
    int glancing_sector_count = 0;
    int direct_sector_count = 0;
    // 0=Direct, 1=Glancing, 2=Blocked, 3=Shielded
    int headline_state = 0;

    // ---- Phase 3A warnings (investor-demo safety rails) ----
    bool warn_fully_blocked = false;
    bool warn_glancing_hold = false;

    double reward = 0.0;
};

class Simulation {
public:
    Simulation();

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;

    // Reactor/Chemistry include internal references; keep Simulation non-movable for determinism.
    Simulation(Simulation&&) = delete;
    Simulation& operator=(Simulation&&) = delete;

    void resetToDataCenterRackScenario();

    // Phase 3A deterministic scenario rig
    void resetToScenario(DemoScenario s);
    void resetToScenario(DemoScenario s, AgentType a);
    void setAgent(AgentType a);
    AgentType agent() const { return agent_type_; }
    void enableCalibrationMode(bool enabled);
    bool calibrationMode() const { return calibration_mode_; }

    // Phase 3B.1 safety harness
    void enableVerificationMode(bool enabled);
    bool verificationMode() const { return verification_mode_; }
    bool runVerificationTest(VerificationTestId id);
    RunSignatures getRunSignatures() const;
    RunSignatures getLastExpectedSignatures() const;
    std::uint32_t getLatestEvents();
    int getTelemetrySamples(TelemetrySampleV1* out_ptr, int cap) const;

    // Phase 3CA: extended monitoring (v2+). v1 remains authoritative for verification.
    void enableTelemetryV2(bool enabled);
    bool telemetryV2Enabled() const { return telemetry_v2_enabled_; }
    int getTelemetrySamplesV2(TelemetrySampleV2* out_ptr, int cap) const;

    // Phase 3CA: autonomy seams (deterministic control input surface)
    void setControlInputs(const ControlInputsV1& in);
    ControlInputsV1 getControlInputs() const { return control_inputs_; }

    // Phase 3CA: deterministic profiling (excluded from verification signatures)
    void enableProfiling(bool enabled);
    bool profilingEnabled() const { return profiling_enabled_; }
    bool getLastProfileSample(ProfileSampleV1* out) const;
    int exportConfigText(char* buf, int cap) const;
    DemoScenario scenario() const { return scenario_; }
    // Scenario-owned nozzle pose: pos in meters, dir must be finite (normalized internally).
    void setNozzlePose(const Vec3d& pos_m, const Vec3d& dir_unit);
    void setNozzleSweepEnabled(bool enabled);

    void commandIgniteOrIncreasePyrolysis();
    void commandStartSuppression();

    // dt must be positive and finite; invalid dt is ignored.
    void step(double dt);

    // Must be side-effect free and NaN-safe for visualizer usage.
    Observation observe() const;

    // Accessor for current nozzle position (meters, world frame)
    const Vec3d& getNozzlePos_m() const { return nozzle_pos_m_; }

    bool isConcluded() const noexcept { return concluded_; }

    // Step 1C: lifecycle state inspection (read-only)
    bool isIgnited() const noexcept { return ignited_; }
    double time_s() const noexcept { return scenario_time_s_; }
    bool isSuppressionEnabled() const { return supp_.config().enabled; }

private:

    static std::vector<Species> buildDefaultSpecies();
    void seedAmbient(Reactor& r);

    ChemistryIndex idx_;

    Reactor reactor_;
    Ventilation vent_;
    Suppression supp_;

    double fuelSolid_kg_ = 50.0;
    double pyrolysis_kgps_ = 0.0;
    double pyrolysisMax_kgps_ = 0.06;

    bool ignited_ = false;
    bool concluded_ = false;

    double lastHRR_W_ = 0.0;
    double inhib_kgm3_ = 0.0;
    double inert_kgm3_ = 0.0;
    double agent_mdot_kgps_ = 0.0;

    // ---- Phase 2A truth telemetry state ----
    double vfep_rpm_ = 0.0;
    double hit_efficiency_0_1_ = 0.0;
    Vec3d  spray_dir_unit_{0.0, 0.0, 0.0};
    Vec3d  draft_vel_mps_{0.0, 0.0, 0.0};
    double jet_momentum_N_ = 0.0;
    double draft_drag_N_   = 0.0;

    // Scenario nozzle axis (unit vector); used as aero module input.
    Vec3d nozzle_dir_unit_{0.0, 0.0, 1.0};

    // ---- Phase 3A scenario rig (sim-owned, deterministic) ----
    DemoScenario scenario_ = DemoScenario::DirectVsGlance;

    // Fire/hotspot truth (meters, world frame)
    Vec3d hotspot_pos_m_{0.0, 0.6, 0.7};

    // Deterministic ignition RNG state (do not use std::rand()).
    std::uint32_t ignition_seed_u32_ = 0u;
    bool ignition_seeded_ = false;

    Vec3d nozzle_pos_m_{-2.0, 1.5, -2.0};
    Vec3d nozzle_dir_unit_scenario_{0.7, -0.15, 0.7};
    bool nozzle_sweep_enabled_ = false;
    double scenario_time_s_ = 0.0;
    double sweep_freq_hz_ = 0.25;
    double sweep_amp_deg_ = 12.0;   

    // Phase 3CA scalability: increase obstacle capacity without changing default behavior.
    // NOTE: legacy scenarios still only populate the first 1–2 obstacles deterministically.
    static constexpr int kMaxObstacles_ = 64;
    AABBd obstacles_[kMaxObstacles_]{};
    int num_obstacles_ = 0;

    // Phase 3CA scalability hook: multi-ray sampling (default 1; verification implicitly uses 1).
    int rays_per_sector_ = 1;

    LiIonRunaway liion_;
    double liionHeat_W_ = 0.0;
    double liionVent_kgps_ = 0.0;

    // Tracks continuous time spent in a "safe" state for UI-friendly termination.
    // Reset to 0 on reset, and whenever the system is outside the safe envelope.
    
    // ---- Phase 2B/2C suppression state (deterministic, time-integrated) ----
    double exposure_kg_ = 0.0;
    double knockdown_target_0_1_ = 0.0;
    double knockdown_0_1_ = 0.0;

    static constexpr int kNumSectors_ = 4;
    std::array<double, kNumSectors_> sector_exposure_kg_{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSectors_> sector_knockdown_target_0_1_{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSectors_> sector_knockdown_0_1_{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSectors_> sector_delivered_mdot_kgps_{{0.0,0.0,0.0,0.0}};

    // ---- Phase 3 geometry suppression state (per sector) ----
    std::array<double, kNumSectors_> sector_occlusion_0_1_{{1.0,1.0,1.0,1.0}};
    std::array<double, kNumSectors_> sector_line_attack_0_1_{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSectors_> sector_net_delivered_mdot_kgps_{{0.0,0.0,0.0,0.0}};

    // ---- Phase 3A stabilization state (per sector) ----
    std::array<double, kNumSectors_> sector_raw_delivered_mdot_kgps_{{0.0,0.0,0.0,0.0}};
    std::array<double, kNumSectors_> sector_shield_0_1_{{1.0,1.0,1.0,1.0}};
    std::array<double, kNumSectors_> occ_stable_0_1_{{1.0,1.0,1.0,1.0}};
    std::array<double, kNumSectors_> shield_stable_0_1_{{1.0,1.0,1.0,1.0}};
    std::array<double, kNumSectors_> loa_smooth_0_1_{{0.0,0.0,0.0,0.0}};
    std::array<int, kNumSectors_> occ_enter_count_{{0,0,0,0}};
    std::array<int, kNumSectors_> occ_exit_count_{{0,0,0,0}};
    std::array<int, kNumSectors_> shield_enter_count_{{0,0,0,0}};
    std::array<int, kNumSectors_> shield_exit_count_{{0,0,0,0}};

    // Tunables
    double loa_power_ = 2.0;
    double loa_min_0_1_ = 0.15;
    double loa_smooth_tau_s_ = 0.20;
    double hysteresis_enter_s_ = 0.10;
    double hysteresis_exit_s_  = 0.20;
    double aabb_pad_m_ = 0.015;

    // Safety rails (time held in condition)
    double fully_blocked_hold_s_ = 0.0;
    double glancing_hold_s_ = 0.0;

    double raw_HRR_W_ = 0.0;
    double effective_HRR_W_ = 0.0;
    SuppressionRegime suppression_regime_ = SuppressionRegime::None;

    

// ---- Phase 3B chemical effectiveness state (deterministic, calibratable) ----
AgentType agent_type_ = AgentType::CleanAgent;
AgentProfile agent_profile_{};

ScenarioFactors scenario_factors_{};

std::array<double, kNumSectors_> sector_utilization_U_0_1_{{0.0,0.0,0.0,0.0}};
std::array<double, kNumSectors_> sector_effective_exposure_kg_{{0.0,0.0,0.0,0.0}};
std::array<double, kNumSectors_> sector_EC50_adj_kg_{{0.0,0.0,0.0,0.0}};
double utilization_U_0_1_ = 0.0;
double effective_exposure_kg_ = 0.0;
double EC50_adj_kg_ = 0.0;

// Calibration mode + signatures
bool calibration_mode_ = false;
double calibration_elapsed_s_ = 0.0;
std::uint32_t run_param_hash_u32_ = 0;
std::uint32_t telemetry_crc_u32_ = 0;

    // ---- Phase 3B.1 verification + telemetry harness ----
    bool verification_mode_ = false;
    RunSignatures run_signatures_{};
    RunSignatures expected_signatures_{};
    std::uint32_t latest_events_bits_ = 0;

    // Telemetry ring buffer (fixed capacity, no dynamic alloc during step)
    static constexpr int kTelemetryCapacity_ = 2048;
    std::array<TelemetrySampleV1, kTelemetryCapacity_> telemetry_rb_{};
    int telemetry_head_ = 0;   // next write
    int telemetry_count_ = 0;  // number valid
    double telemetry_dt_s_ = 0.10;
    double telemetry_next_t_s_ = 0.0;

    // ---- Phase 3CA: Telemetry v2 ring buffer (derived monitoring; never hashed) ----
    bool telemetry_v2_enabled_ = false;
    static constexpr int kTelemetryV2Capacity_ = 2048;
    std::array<TelemetrySampleV2, kTelemetryV2Capacity_> telemetry_v2_rb_{};
    int telemetry_v2_head_ = 0;
    int telemetry_v2_count_ = 0;

    // ---- Phase 3CA: autonomy seams ----
    ControlInputsV1 control_inputs_{};

    // ---- Phase 3CA: profiling (excluded from verification hashes) ----
    bool profiling_enabled_ = false;
    ProfileSampleV1 last_profile_{};

    // Event edge detection state
    bool prev_occluded_any_ = false;
    bool prev_kd_ge_0p5_ = false;
    bool prev_kd_ge_0p9_ = false;
    bool prev_hrr_below_100kW_ = false;
    double prev_exposure_kg_ = 0.0;
    double prev_effective_exposure_kg_ = 0.0;
    double prev_kd_target_0_1_ = 0.0;

    // Tunables (simple, monotonic, deterministic)
    double exposure_half_kg_ = 0.35;
    double exposure_hill_n_  = 2.0;
    double tau_knockdown_rise_s_ = 2.0;
    double tau_knockdown_fall_s_ = 4.0;
    bool   enable_recovery_ = true;
    double tau_exposure_decay_s_ = 10.0;

    double safeHold_s_ = 0.0;
};

} // namespace vfep
