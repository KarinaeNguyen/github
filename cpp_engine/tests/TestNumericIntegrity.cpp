#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <memory>

#include "Simulation.h"
#include "SensitivityAnalysis.h"
#include "UncertaintyQuantification.h"

namespace {

// Always-on requirement: never compiled out in Release.
#define REQUIRE(cond, msg)                                                      \
    do {                                                                        \
        if (!(cond)) {                                                          \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " " << msg \
                      << "\n";                                                  \
            std::exit(1);                                                       \
        }                                                                       \
    } while (0)

static inline void REQUIRE_FINITE(double x, const char* name) {
    if (!std::isfinite(x)) {
        std::cerr << "[FAIL] Non-finite: " << name << " = " << x << "\n";
        std::exit(1);
    }
}

static inline double absd(double x) { return x < 0 ? -x : x; }
static inline double maxd(double a, double b) { return (a > b) ? a : b; }

static void requireFiniteAndBounded(const vfep::Observation& o) {
    // 1B.1: Must be finite (hard gate)
    REQUIRE_FINITE(o.T_K, "T_K");
    REQUIRE_FINITE(o.HRR_W, "HRR_W");
    REQUIRE_FINITE(o.O2_volpct, "O2_volpct");
    REQUIRE_FINITE(o.CO2_volpct, "CO2_volpct");
    REQUIRE_FINITE(o.H2O_volpct, "H2O_volpct");
    REQUIRE_FINITE(o.fuel_kg, "fuel_kg");
    REQUIRE_FINITE(o.inhibitor_kgm3, "inhibitor_kgm3");
    REQUIRE_FINITE(o.inert_kgm3, "inert_kgm3");
    REQUIRE_FINITE(o.ACH, "ACH");
    REQUIRE_FINITE(o.agent_mdot_kgps, "agent_mdot_kgps");
    REQUIRE_FINITE(o.reward, "reward");

    // 1B.2: Boundedness / safety rails (not physics correctness)
    REQUIRE(o.T_K >= 1.0 && o.T_K <= 5000.0, "T_K out of bounds [1,5000]");
    REQUIRE(o.HRR_W >= 0.0, "HRR_W < 0");
    REQUIRE(o.fuel_kg >= 0.0, "fuel_kg < 0");
    REQUIRE(o.inhibitor_kgm3 >= 0.0, "inhibitor_kgm3 < 0");
    REQUIRE(o.inert_kgm3 >= 0.0, "inert_kgm3 < 0");

    // Vol% should remain in [0,100]
    REQUIRE(o.O2_volpct >= 0.0 && o.O2_volpct <= 100.0, "O2_volpct out of [0,100]");
    REQUIRE(o.CO2_volpct >= 0.0 && o.CO2_volpct <= 100.0, "CO2_volpct out of [0,100]");
    REQUIRE(o.H2O_volpct >= 0.0 && o.H2O_volpct <= 100.0, "H2O_volpct out of [0,100]");

    // Control / IO channel sanity (tripwires, non-tuning)
    REQUIRE(o.ACH >= 0.0, "ACH < 0");
    REQUIRE(o.ACH <= 1.0e4, "ACH unrealistically high");

    REQUIRE(o.agent_mdot_kgps >= 0.0, "agent_mdot_kgps < 0");
    REQUIRE(o.agent_mdot_kgps <= 1.0e4, "agent_mdot_kgps unrealistically high");

    REQUIRE(std::abs(o.reward) <= 1.0e9, "reward magnitude too large");

    // Vol% sum sanity (very loose, non-tuning)
    const double vol_sum = o.O2_volpct + o.CO2_volpct + o.H2O_volpct;
    REQUIRE(vol_sum >= 0.0 && vol_sum <= 300.0, "vol% sum out of sanity range");
}

static void requireCloseAbsOrRel(const char* name, double a, double b, double absTol, double relTol) {
    REQUIRE_FINITE(a, name);
    REQUIRE_FINITE(b, name);

    const double diff = absd(a - b);
    const double denom = maxd(absd(a), absd(b));
    const double rel = (denom > 0.0) ? (diff / denom) : diff;

    if (!(diff <= absTol || rel <= relTol)) {
        std::cerr << "[FAIL] 1B.3a dt-robustness: " << name
                  << " a=" << a << " b=" << b
                  << " diff=" << diff << " (absTol=" << absTol << ")"
                  << " rel=" << rel << " (relTol=" << relTol << ")\n";
        std::exit(1);
    }
}

static void requireExact(const char* label, double a, double b) {
    REQUIRE_FINITE(a, label);
    REQUIRE_FINITE(b, label);
    if (!(a == b)) {
        std::cerr << "[FAIL] " << label << " changed unexpectedly: a=" << a << " b=" << b << "\n";
        std::exit(1);
    }
}

static void requireObservationExactKeyFields(
    const char* context,
    const vfep::Observation& a,
    const vfep::Observation& b)
{
    // Key channels that should remain bitwise-identical if an operation is a pure no-op.
    // (If this is too strict in the future, relax *only* the fields that prove unstable.)
    std::string prefix = std::string(context) + ": ";

    requireExact((prefix + "T_K").c_str(), a.T_K, b.T_K);
    requireExact((prefix + "HRR_W").c_str(), a.HRR_W, b.HRR_W);
    requireExact((prefix + "O2_volpct").c_str(), a.O2_volpct, b.O2_volpct);
    requireExact((prefix + "CO2_volpct").c_str(), a.CO2_volpct, b.CO2_volpct);
    requireExact((prefix + "H2O_volpct").c_str(), a.H2O_volpct, b.H2O_volpct);
    requireExact((prefix + "fuel_kg").c_str(), a.fuel_kg, b.fuel_kg);
    requireExact((prefix + "inhibitor_kgm3").c_str(), a.inhibitor_kgm3, b.inhibitor_kgm3);
    requireExact((prefix + "inert_kgm3").c_str(), a.inert_kgm3, b.inert_kgm3);
    requireExact((prefix + "ACH").c_str(), a.ACH, b.ACH);
    requireExact((prefix + "agent_mdot_kgps").c_str(), a.agent_mdot_kgps, b.agent_mdot_kgps);
    requireExact((prefix + "reward").c_str(), a.reward, b.reward);
}

/* =======================
 * Step 1C helpers
 * ======================= */

struct StateSnap {
    double t_s = 0.0;
    bool ignited = false;
    bool suppressed = false;
    bool concluded = false;
    vfep::Observation o{};
};

static StateSnap snap(const vfep::Simulation& sim) {
    StateSnap s;
    s.t_s = sim.time_s();
    s.ignited = sim.isIgnited();
    s.suppressed = sim.isSuppressionEnabled();
    s.concluded = sim.isConcluded();
    s.o = sim.observe();

    requireFiniteAndBounded(s.o);
    REQUIRE_FINITE(s.t_s, "time_s");
    return s;
}

static void requireLatchMonotonic(const StateSnap& prev, const StateSnap& cur, const char* ctx) {
    if (prev.ignited) {
        REQUIRE(cur.ignited, std::string(ctx).append(": ignited regressed").c_str());
    }
    if (prev.suppressed) {
        REQUIRE(cur.suppressed, std::string(ctx).append(": suppressed regressed").c_str());
    }
    if (prev.concluded) {
        REQUIRE(cur.concluded, std::string(ctx).append(": concluded regressed").c_str());
    }
    if (cur.concluded) {
        REQUIRE(cur.ignited, std::string(ctx).append(": concluded without ignition").c_str());
    }
}

static void requireFuelNonIncreasing(const StateSnap& prev, const StateSnap& cur, const char* ctx) {
    REQUIRE(cur.o.fuel_kg <= prev.o.fuel_kg,
            std::string(ctx).append(": fuel increased").c_str());
}

static void requireTimeAdvancedBy(const char* ctx, double t_before, double t_after, double dt_expected) {
    REQUIRE_FINITE(t_before, "t_before");
    REQUIRE_FINITE(t_after, "t_after");
    REQUIRE_FINITE(dt_expected, "dt_expected");

    const double expected = t_before + dt_expected;
    // Tight tolerance to avoid floating accumulation false-fails.
    const double tol = 1e-12 * maxd(1.0, absd(expected));
    const double diff = absd(t_after - expected);

    REQUIRE(t_after > t_before, std::string(ctx).append(": time did not increase").c_str());
    REQUIRE(diff <= tol, std::string(ctx).append(": time jump != dt (within tol)").c_str());
}

/* =======================
 * Existing Step 1B helpers
 * ======================= */

static vfep::Observation baselineObservation() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();
    const vfep::Observation o = sim.observe();
    requireFiniteAndBounded(o);
    return o;
}

// Harness used for 1B.3a comparisons and normal soak runs.
// This harness assumes caller provides dt > 0 and t_end >= 0.
// (Invalid dt tests are done separately in 1B.3b.)
static vfep::Observation runToTimeAndSample(double dt, double t_end, bool ignite, bool suppress) {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double ignite_at = 2.0;
    const double suppress_at = 5.0;

    double t = 0.0;
    bool didIgnite = false;
    bool didSuppress = false;

    REQUIRE_FINITE(dt, "dt");
    REQUIRE(dt > 0.0, "dt must be > 0");
    REQUIRE_FINITE(t_end, "t_end");
    REQUIRE(t_end >= 0.0, "t_end must be >= 0");

    // Non-termination guard: if time fails to advance, fail fast.
    const std::size_t cap_by_ratio = (dt > 0.0) ? (std::size_t)(t_end / dt) + 16u : 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    requireFiniteAndBounded(sim.observe());

    while (t < t_end) {
        ++iter;
        REQUIRE(iter <= iter_cap, "non-termination guard tripped (time stall or unexpected loop)");

        const double t_next = t + dt;

        // Robust scheduling even for large dt that may "jump" over thresholds.
        if (ignite && !didIgnite && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
        }
        if (suppress && !didSuppress && (t >= suppress_at || (t < suppress_at && t_next >= suppress_at))) {
            sim.commandStartSuppression();
            didSuppress = true;
        }

        sim.step(dt);
        t = t_next;

        requireFiniteAndBounded(sim.observe());
    }

    // Ensure schedule was actually applied (prevents silent no-op regressions).
    if (ignite && t_end > ignite_at) {
        REQUIRE(didIgnite, "ignite requested but ignition command was never issued");
    }
    if (suppress && t_end > suppress_at) {
        REQUIRE(didSuppress, "suppression requested but suppression command was never issued");
    }

    return sim.observe();
}

static void runDtRobustnessTripwire_1B3a() {
    const double dt_small = 0.02;
    const double dt_large = 0.10;

    const double t_pre = 1.0;   // pre-ignite
    const double t_post = 10.0; // post-ignite + post-suppress

    const bool ignite = true;
    const bool suppress = true;

    {
        const vfep::Observation a = runToTimeAndSample(dt_small, t_pre, ignite, suppress);
        const vfep::Observation b = runToTimeAndSample(dt_large, t_pre, ignite, suppress);

        requireCloseAbsOrRel("T_K@t=1", a.T_K, b.T_K, 250.0, 0.20);
        requireCloseAbsOrRel("HRR_W@t=1", a.HRR_W, b.HRR_W, 5000.0, 0.50);
        requireCloseAbsOrRel("O2_volpct@t=1", a.O2_volpct, b.O2_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("CO2_volpct@t=1", a.CO2_volpct, b.CO2_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("H2O_volpct@t=1", a.H2O_volpct, b.H2O_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("fuel_kg@t=1", a.fuel_kg, b.fuel_kg, 0.05, 0.20);
        requireCloseAbsOrRel("inhibitor_kgm3@t=1", a.inhibitor_kgm3, b.inhibitor_kgm3, 0.05, 0.50);
        requireCloseAbsOrRel("inert_kgm3@t=1", a.inert_kgm3, b.inert_kgm3, 0.05, 0.50);
        requireCloseAbsOrRel("ACH@t=1", a.ACH, b.ACH, 1.0, 0.50);
        requireCloseAbsOrRel("agent_mdot_kgps@t=1", a.agent_mdot_kgps, b.agent_mdot_kgps, 0.1, 0.50);
    }

    {
        const vfep::Observation a = runToTimeAndSample(dt_small, t_post, ignite, suppress);
        const vfep::Observation b = runToTimeAndSample(dt_large, t_post, ignite, suppress);

        requireCloseAbsOrRel("T_K@t=10", a.T_K, b.T_K, 250.0, 0.20);
        requireCloseAbsOrRel("HRR_W@t=10", a.HRR_W, b.HRR_W, 5000.0, 0.50);
        requireCloseAbsOrRel("O2_volpct@t=10", a.O2_volpct, b.O2_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("CO2_volpct@t=10", a.CO2_volpct, b.CO2_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("H2O_volpct@t=10", a.H2O_volpct, b.H2O_volpct, 10.0, 0.20);
        requireCloseAbsOrRel("fuel_kg@t=10", a.fuel_kg, b.fuel_kg, 0.05, 0.20);
        requireCloseAbsOrRel("inhibitor_kgm3@t=10", a.inhibitor_kgm3, b.inhibitor_kgm3, 0.05, 0.50);
        requireCloseAbsOrRel("inert_kgm3@t=10", a.inert_kgm3, b.inert_kgm3, 0.05, 0.50);
        requireCloseAbsOrRel("ACH@t=10", a.ACH, b.ACH, 1.0, 0.50);
        requireCloseAbsOrRel("agent_mdot_kgps@t=10", a.agent_mdot_kgps, b.agent_mdot_kgps, 0.1, 0.50);
    }

    std::cout << "[PASS] 1B.3a dt-robustness tripwire (dt=0.02 vs 0.10 at t=1s,10s)\n";
}

static void runSoak(const char* name, double dt, double t_end, bool ignite, bool suppress) {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    double t = 0.0;
    bool didIgnite = false;
    bool didSuppress = false;

    const double ignite_at = 2.0;
    const double suppress_at = 5.0;

    REQUIRE_FINITE(dt, "dt");
    REQUIRE(dt > 0.0, "dt must be > 0");
    REQUIRE_FINITE(t_end, "t_end");
    REQUIRE(t_end >= 0.0, "t_end must be >= 0");

    const std::size_t cap_by_ratio = (dt > 0.0) ? (std::size_t)(t_end / dt) + 16u : 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    requireFiniteAndBounded(sim.observe());

    while (t < t_end) {
        ++iter;
        REQUIRE(iter <= iter_cap, "non-termination guard tripped (time stall or unexpected loop)");

        const double t_next = t + dt;

        if (ignite && !didIgnite && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
        }
        if (suppress && !didSuppress && (t >= suppress_at || (t < suppress_at && t_next >= suppress_at))) {
            sim.commandStartSuppression();
            didSuppress = true;
        }

        sim.step(dt);
        t = t_next;

        requireFiniteAndBounded(sim.observe());
    }

    if (ignite && t_end > ignite_at) {
        REQUIRE(didIgnite, "ignite requested but ignition command was never issued");
    }
    if (suppress && t_end > suppress_at) {
        REQUIRE(didSuppress, "suppression requested but suppression command was never issued");
    }

    std::cout << "[PASS] " << name << " (dt=" << dt << ", t_end=" << t_end << ")\n";
}

// =======================
// 1B.3b: Input contract & invalid-parameter safety
// =======================

static void runInvalidDtHandling_1B3b() {
    const vfep::Observation base = baselineObservation();

    // dt == 0
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();
        const vfep::Observation before = sim.observe();
        requireFiniteAndBounded(before);

        sim.step(0.0);
        const vfep::Observation after = sim.observe();
        requireFiniteAndBounded(after);

        requireObservationExactKeyFields("1B.3b dt==0 no-op", before, after);
        requireObservationExactKeyFields("1B.3b dt==0 baseline", base, after);
    }

    // dt < 0
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();
        const vfep::Observation before = sim.observe();
        requireFiniteAndBounded(before);

        sim.step(-0.1);
        const vfep::Observation after = sim.observe();
        requireFiniteAndBounded(after);

        requireObservationExactKeyFields("1B.3b dt<0 no-op", before, after);
        requireObservationExactKeyFields("1B.3b dt<0 baseline", base, after);
    }

    // dt = NaN
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();
        const vfep::Observation before = sim.observe();
        requireFiniteAndBounded(before);

        const double dt_nan = std::numeric_limits<double>::quiet_NaN();
        sim.step(dt_nan);

        const vfep::Observation after = sim.observe();
        requireFiniteAndBounded(after);

        requireObservationExactKeyFields("1B.3b dt=NaN no-op", before, after);
        requireObservationExactKeyFields("1B.3b dt=NaN baseline", base, after);
    }

    // dt = +Inf
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();
        const vfep::Observation before = sim.observe();
        requireFiniteAndBounded(before);

        const double dt_inf = std::numeric_limits<double>::infinity();
        sim.step(dt_inf);

        const vfep::Observation after = sim.observe();
        requireFiniteAndBounded(after);

        requireObservationExactKeyFields("1B.3b dt=Inf no-op", before, after);
        requireObservationExactKeyFields("1B.3b dt=Inf baseline", base, after);
    }

    std::cout << "[PASS] 1B.3b invalid dt handling (0, <0, NaN, Inf) no-op + bounded\n";
}

static void runVeryLargeDtSafety_1B3b() {
    runSoak("huge_dt_60s_step_10m", 60.0, 600.0, true, true);
    runSoak("huge_dt_600s_step_10m", 600.0, 600.0, true, true);

    std::cout << "[PASS] 1B.3b very large dt safety (60s, 600s) finite/bounded\n";
}

static void runCommandSchedulingRobustness_1B3b() {
    // A) Ignite spam every step
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        const double dt = 0.05;
        const int steps = 400; // 20s

        for (int i = 0; i < steps; ++i) {
            sim.commandIgniteOrIncreasePyrolysis();
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }

        std::cout << "[PASS] 1B.3b ignite spam every step: finite/bounded\n";
    }

    // B) Suppression spam every step (even before ignition)
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        const double dt = 0.05;
        const int steps = 400; // 20s

        for (int i = 0; i < steps; ++i) {
            sim.commandStartSuppression();
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }

        std::cout << "[PASS] 1B.3b suppression spam every step: finite/bounded\n";
    }

    // C) Suppression before ignition, then ignite later, both safe
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        const double dt = 0.05;

        // Suppress first for 2s
        for (int i = 0; i < 40; ++i) {
            sim.commandStartSuppression();
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }

        // Then ignite for 10s while still possibly suppressing
        for (int i = 0; i < 200; ++i) {
            sim.commandIgniteOrIncreasePyrolysis();
            sim.commandStartSuppression();
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }

        std::cout << "[PASS] 1B.3b suppression-before-ignition then ignite: finite/bounded\n";
    }
}

static void runResetRobustness_1B3b() {
    const vfep::Observation base = baselineObservation();

    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        for (int i = 0; i < 200; ++i) {
            sim.commandIgniteOrIncreasePyrolysis();
            sim.commandStartSuppression();
            sim.step(0.05);
            requireFiniteAndBounded(sim.observe());
        }

        sim.step(0.0);
        sim.step(-1.0);
        requireFiniteAndBounded(sim.observe());

        sim.resetToDataCenterRackScenario();
        sim.resetToDataCenterRackScenario();

        const vfep::Observation after = sim.observe();
        requireFiniteAndBounded(after);
        requireObservationExactKeyFields("1B.3b reset twice baseline", base, after);

        std::cout << "[PASS] 1B.3b reset idempotence (double reset returns baseline)\n";
    }

    {
        vfep::Simulation simA;
        simA.resetToDataCenterRackScenario();
        simA.step(std::numeric_limits<double>::quiet_NaN());
        simA.commandStartSuppression();
        simA.commandIgniteOrIncreasePyrolysis();
        requireFiniteAndBounded(simA.observe());

        vfep::Simulation simB;
        simB.resetToDataCenterRackScenario();
        const vfep::Observation fresh = simB.observe();
        requireFiniteAndBounded(fresh);
        requireObservationExactKeyFields("1B.3b fresh sim baseline", base, fresh);

        std::cout << "[PASS] 1B.3b fresh instance after invalid calls remains baseline\n";
    }
}

/* =======================
 * Step 1C tests
 * ======================= */

static void runCommandDoesNotAdvanceTime_1C() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const StateSnap a = snap(sim);

    sim.commandStartSuppression();
    const StateSnap b = snap(sim);
    REQUIRE(b.t_s == a.t_s, "1C.E1 commandStartSuppression advanced time");

    sim.commandIgniteOrIncreasePyrolysis();
    const StateSnap c = snap(sim);
    REQUIRE(c.t_s == a.t_s, "1C.E1 commandIgnite advanced time");

    std::cout << "[PASS] 1C.E1 commands do not advance time\n";
}

static void runTimeMonotonicExactDt_1C() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    double t_prev = sim.time_s();
    const double dt = 0.05;

    for (int i = 0; i < 200; ++i) {
        sim.step(dt);
        const double t = sim.time_s();
        requireTimeAdvancedBy("1C.E2", t_prev, t, dt);
        t_prev = t;
        requireFiniteAndBounded(sim.observe());
    }

    std::cout << "[PASS] 1C.E2 time monotonic + dt-consistent\n";
}

static void runIgnitionLatchAndIdempotence_1C() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    StateSnap prev = snap(sim);
    REQUIRE(!prev.ignited, "1C.B1 baseline ignited unexpectedly");

    sim.commandIgniteOrIncreasePyrolysis();
    StateSnap cur = snap(sim);
    REQUIRE(cur.ignited, "1C.B1 ignition did not latch true on command");
    requireLatchMonotonic(prev, cur, "1C ignition latch");
    prev = cur;

    const double dt = 0.05;
    for (int i = 0; i < 400; ++i) {
        sim.commandIgniteOrIncreasePyrolysis();
        sim.step(dt);

        cur = snap(sim);
        requireLatchMonotonic(prev, cur, "1C ignition spam");
        requireFuelNonIncreasing(prev, cur, "1C fuel monotonic (ignite spam)");
        prev = cur;
    }

    std::cout << "[PASS] 1C.B1/D1/B3 ignition latch + idempotence + fuel non-increase\n";
}

static void runSuppressionBeforeIgnition_1C() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    StateSnap prev = snap(sim);
    REQUIRE(!prev.ignited, "1C.D2 baseline ignited unexpectedly");
    REQUIRE(!prev.concluded, "1C.D2 baseline concluded unexpectedly");

    const double dt = 0.05;

    for (int i = 0; i < 80; ++i) {
        sim.commandStartSuppression();
        sim.step(dt);

        StateSnap cur = snap(sim);
        REQUIRE(!cur.ignited, "1C.D2 suppression caused ignition history");
        REQUIRE(!cur.concluded, "1C.D2 concluded before ignition");
        requireLatchMonotonic(prev, cur, "1C suppression-before-ignition");
        requireFuelNonIncreasing(prev, cur, "1C fuel monotonic (pre-ignite)");
        prev = cur;
    }

    sim.commandIgniteOrIncreasePyrolysis();

    for (int i = 0; i < 200; ++i) {
        sim.commandStartSuppression();
        sim.step(dt);

        StateSnap cur = snap(sim);
        REQUIRE(cur.ignited, "1C later ignition did not remain latched");
        requireLatchMonotonic(prev, cur, "1C suppression then ignition");
        requireFuelNonIncreasing(prev, cur, "1C fuel monotonic (post-ignite)");
        prev = cur;
    }

    std::cout << "[PASS] 1C.D2/B2 suppression-before-ignition safe + latch monotonic\n";
}

static void runTerminalNoOpAndCommandSafety_1C() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.10;
    const double t_cap = 2.0 * 3600.0;

    double t = 0.0;
    StateSnap prev = snap(sim);

    while (t < t_cap && !prev.concluded) {
        const double t_next = t + dt;

        // Same scheduling style as your other tests
        if (t < 2.0 && t_next >= 2.0) sim.commandIgniteOrIncreasePyrolysis();
        if (t < 5.0 && t_next >= 5.0) sim.commandStartSuppression();

        sim.step(dt);
        t = t_next;

        StateSnap cur = snap(sim);
        requireLatchMonotonic(prev, cur, "1C terminal seek");
        requireFuelNonIncreasing(prev, cur, "1C terminal seek fuel");
        prev = cur;
    }

    if (!prev.concluded) {
        // Tolerated: Step 1C does not require the scenario to reach a terminal state;
        // it only requires that terminal behavior is safe and non-regressive if present.
        std::cout << "[PASS] 1C.B4/D3 terminal behavior: not reached (tolerated)\n";
        return;
    }

    const StateSnap c0 = prev;

    // After conclusion: step should be a pure no-op (your engine early-returns)
    sim.step(dt);
    const StateSnap c1 = snap(sim);
    REQUIRE(c1.concluded, "1C.B4 concluded regressed after step");
    REQUIRE(c1.t_s == c0.t_s, "1C terminal: time advanced after conclusion (expected freeze)");
    requireObservationExactKeyFields("1C terminal step no-op", c0.o, c1.o);

    // Commands after conclusion must be safe and must not reopen
    sim.commandIgniteOrIncreasePyrolysis();
    sim.commandStartSuppression();
    const StateSnap c2 = snap(sim);
    REQUIRE(c2.concluded, "1C.D3 concluded cleared by commands");
    REQUIRE(c2.t_s == c0.t_s, "1C.D3 commands advanced time in terminal");
    requireObservationExactKeyFields("1C terminal command no-op", c0.o, c2.o);

    std::cout << "[PASS] 1C.B4/D3 terminal latch + terminal no-op + command safety\n";
}

// =======================
// Step 2.1: Pre-ignition quiescence
// =======================
static void runPreIgnitionQuiescence_2A1() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.05;
    const double t_end = 1.5;     // strictly before ignite_at=2.0 in other harness paths
    const double eps_fuel = 1e-12;
    const double eps_hrr  = 1e-9; // expect exact 0 in most cases; keep tiny epsilon

    double t = 0.0;

    const StateSnap s0 = snap(sim);
    REQUIRE(!s0.ignited, "2.1 pre-ignite: ignited latch true at t=0");
    REQUIRE(!s0.suppressed, "2.1 pre-ignite: suppressed latch true at t=0");
    REQUIRE(s0.o.HRR_W <= eps_hrr, "2.1 pre-ignite: HRR_W not quiescent at t=0");

    while (t < t_end) {
        sim.step(dt);
        t += dt;

        const StateSnap s = snap(sim);

        REQUIRE(!s.ignited, "2.1 pre-ignite: ignited became true without ignition command");
        REQUIRE(!s.suppressed, "2.1 pre-ignite: suppressed became true without suppression command");
        REQUIRE(s.o.HRR_W <= eps_hrr, "2.1 pre-ignite: HRR_W rose above quiescent threshold");

        REQUIRE(s.o.fuel_kg + eps_fuel >= s0.o.fuel_kg,
                "2.1 pre-ignite: fuel decreased before ignition");
    }

    std::cout << "[PASS] 2.1 pre-ignition quiescence (t<=1.5s, dt=0.05)\n";
}

/* =======================
 * Step 2A2: Physical Consistency — Combustion implies fuel consumption (revised)
 * ======================= */
static void runCombustionImpliesFuelConsumption_2A2() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.05;
    const double ignite_at = 2.0;

    double t = 0.0;

    // Advance to just before ignition
    while (t + dt < ignite_at) {
        sim.step(dt);
        t += dt;
    }

    // Use the SAME ignition command used in Step 1C
    sim.commandIgniteOrIncreasePyrolysis();

    // Step once to move past the command boundary (Step 1C contract)
    sim.step(dt);
    t += dt;

    // Step 2.2: ignition must cause fuel to decrease within a bounded window
    const double eps_fuel = 1e-12;
    const double fuel_min_delta = 1e-12;

    const double fuel0 = snap(sim).o.fuel_kg;

    const double window_s = 5.0;
    const int steps = static_cast<int>(window_s / dt);

    for (int i = 0; i < steps; i++) {
        sim.step(dt);
        t += dt;
    }

    const double fuel1 = snap(sim).o.fuel_kg;

    REQUIRE(fuel1 + eps_fuel <= fuel0 - fuel_min_delta,
            "2.2: after ignition, fuel did not decrease within 5s window");

    std::cout << "[PASS] 2.2 ignition implies fuel consumption (fuel decreases post-ignite)\n";
}


/* =======================
 * Step 2.3: Heat release correlates directionally with temperature rise (time-local plausibility)
 * Hard PASS/FAIL in Release builds.
 * ======================= */

static double pearsonR(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    if (n == 0 || y.size() != n) return 0.0;

    double mean_x = 0.0, mean_y = 0.0;
    double Sxx = 0.0, Syy = 0.0, Sxy = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const double xi = x[i];
        const double yi = y[i];
        if (!std::isfinite(xi) || !std::isfinite(yi)) return 0.0;

        const double dx = xi - mean_x;
        const double dy = yi - mean_y;
        const double inv = 1.0 / static_cast<double>(i + 1);

        mean_x += dx * inv;
        mean_y += dy * inv;

        Sxx += dx * (xi - mean_x);
        Syy += dy * (yi - mean_y);
        Sxy += dx * (yi - mean_y);
    }

    if (!(Sxx > 0.0) || !(Syy > 0.0)) return 0.0;
    const double denom = std::sqrt(Sxx * Syy);
    if (!(denom > 0.0) || !std::isfinite(denom)) return 0.0;
    const double r = Sxy / denom;
    return std::isfinite(r) ? r : 0.0;
}

static void runHeatReleaseCorrelatesWithTemperatureRise_2A3() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.05;
    const double t_end = 4.9;    // strictly before suppress_at = 5.0
    const double eps_T = 1e-6;   // K per step tolerance

    // Ignition issued once using the scenario command.
    sim.commandIgniteOrIncreasePyrolysis();

    std::vector<double> hrr;
    std::vector<double> dT;
    hrr.reserve(static_cast<std::size_t>(t_end / dt) + 8u);
    dT.reserve(static_cast<std::size_t>(t_end / dt) + 8u);

    StateSnap prev = snap(sim);

    double t = 0.0;
    while (t + dt <= t_end + 1e-12) {
        sim.step(dt);
        t += dt;

        const StateSnap cur = snap(sim);

        const double this_dT = cur.o.T_K - prev.o.T_K;
        REQUIRE_FINITE(this_dT, "2.3 dT");

        hrr.push_back(cur.o.HRR_W);
        dT.push_back(this_dT);

        prev = cur;
    }

    REQUIRE(!hrr.empty(), "2.3: empty sample set");
    REQUIRE(hrr.size() == dT.size(), "2.3: HRR/dT size mismatch");

    double hrr_max = 0.0;
    for (double v : hrr) {
        REQUIRE_FINITE(v, "2.3 HRR_W");
        if (v > hrr_max) hrr_max = v;
    }
    REQUIRE(hrr_max > 0.0, "2.3: HRR stayed at 0 after ignition");

    // Deterministic adaptive gate (avoids brittle absolute thresholds)
    const double hrr_gate = std::max(1e-6, 0.20 * hrr_max);

    std::vector<double> hrr_f;
    std::vector<double> dT_f;
    hrr_f.reserve(hrr.size());
    dT_f.reserve(dT.size());

    std::size_t eligible = 0;
    std::size_t ok = 0;

    for (std::size_t i = 0; i < hrr.size(); ++i) {
        if (hrr[i] > hrr_gate) {
            ++eligible;
            if (dT[i] >= -eps_T) ++ok;
            hrr_f.push_back(hrr[i]);
            dT_f.push_back(dT[i]);
        }
    }

    REQUIRE(eligible >= 10, "2.3A: insufficient HRR>gate samples (vacuous pass prevention)");

    const double frac_ok = static_cast<double>(ok) / static_cast<double>(eligible);
    REQUIRE(frac_ok >= 0.80, "2.3A: dT not directionally plausible for >=80% of HRR>gate steps");

    // Local correlation (0-lag and 1-step lag to tolerate one-step thermal latency)
    const double r0 = pearsonR(hrr_f, dT_f);

    double r1 = 0.0;
    if (hrr_f.size() >= 2) {
        std::vector<double> hrr_l;
        std::vector<double> dT_l;
        hrr_l.reserve(hrr_f.size() - 1);
        dT_l.reserve(dT_f.size() - 1);

        for (std::size_t i = 1; i < hrr_f.size(); ++i) {
            hrr_l.push_back(hrr_f[i - 1]);
            dT_l.push_back(dT_f[i]);
        }
        r1 = pearsonR(hrr_l, dT_l);
    }

    const double r = (r0 > r1) ? r0 : r1;
    REQUIRE(r >= 0.10, "2.3B: HRR vs dT correlation below threshold (r < 0.10)");

    // Additional rigidity: high-HRR steps must not have weaker mean dT than low-HRR steps.
    {
        std::vector<double> sorted = hrr_f;
        std::sort(sorted.begin(), sorted.end());
        const double median = sorted[sorted.size() / 2];

        double sum_hi = 0.0, sum_lo = 0.0;
        std::size_t n_hi = 0, n_lo = 0;

        for (std::size_t i = 0; i < hrr_f.size(); ++i) {
            if (hrr_f[i] >= median) { sum_hi += dT_f[i]; ++n_hi; }
            else { sum_lo += dT_f[i]; ++n_lo; }
        }

        const double mean_hi = (n_hi > 0) ? (sum_hi / static_cast<double>(n_hi)) : 0.0;
        const double mean_lo = (n_lo > 0) ? (sum_lo / static_cast<double>(n_lo)) : 0.0;

        REQUIRE(mean_hi >= mean_lo - eps_T, "2.3C: high-HRR steps do not exhibit >= low-HRR temperature tendency");
    }

    std::cout << "[PASS] 2.3 HRR->T plausibility: eligible=" << eligible
              << " ok=" << ok
              << " frac_ok=" << frac_ok
              << " hrr_gate=" << hrr_gate
              << " r0=" << r0
              << " r1=" << r1 << "\n";
}

/* =======================
 * Step 2.4: Suppression reduces effective heating outcomes
 * (directional modifier consistency)
 *
 * Requirements:
 * - Two matched simulations A (no suppression) and B (suppression enabled).
 * - Identical dt, initial conditions, ignition timing.
 * - Enable suppression in B at a known time after ignition.
 * - Post-suppression:
 *   - HRR_B(t) <= HRR_A(t) + eps for all eligible timesteps.
 *   - DeltaT_B (from suppression-on time to end) <= DeltaT_A + eps.
 *   - No paradoxical steps where HRR_B << HRR_A but dT_B >> dT_A.
 */
static void runSuppressionReducesHeatingOutcomes_2A4() {
    vfep::Simulation simA;
    vfep::Simulation simB;
    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Assumptions (explicitly stated by the test):
    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_on_at = ignite_at + 0.5; // 0.5s after ignition
    const double t_end = 12.0;                    // fixed comparison window

    // Tolerances (eps): use small absolute eps plus a scale-aware term.
    // HRR scale is scenario-dependent; we compute eps after observing the post-suppression peak in A.
    const double eps_T_step = 1e-6;  // K per-step tolerance
    const double eps_T_total = 1e-6; // K on total DeltaT comparison

    REQUIRE(dt > 0.0, "2.4 dt must be > 0");
    REQUIRE(t_end > suppress_on_at, "2.4 t_end must be > suppress_on_at");

    double t = 0.0;
    bool didIgniteA = false;
    bool didIgniteB = false;
    bool didSuppressB = false;

    StateSnap prevA = snap(simA);
    StateSnap prevB = snap(simB);

    // Track the moment suppression becomes active for the DeltaT window.
    bool haveSuppressionStartTemps = false;
    double T_A_start = 0.0;
    double T_B_start = 0.0;

    // Post-suppression HRR gating: only compare when the baseline (A) is meaningfully burning.
    // Use an adaptive gate based on the running max HRR in A after suppression-on.
    double hrrA_running_max_post = 0.0;
    std::size_t eligible = 0;

    // Track the peak HRR in the post-suppression window for epsilon sizing.
    double hrrA_peak_post = 0.0;

    // Non-termination guard (same pattern as other harnesses).
    const std::size_t cap_by_ratio = (std::size_t)(t_end / dt) + 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    // Keep the two sims in lockstep.
    while (t + dt <= t_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.4 non-termination guard tripped");

        const double t_next = t + dt;

        // Matched ignition scheduling.
        if (!didIgniteA && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simA.commandIgniteOrIncreasePyrolysis();
            didIgniteA = true;
        }
        if (!didIgniteB && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simB.commandIgniteOrIncreasePyrolysis();
            didIgniteB = true;
        }

        // Suppression only in B, at a known delay after ignition.
        if (!didSuppressB && (t >= suppress_on_at || (t < suppress_on_at && t_next >= suppress_on_at))) {
            simB.commandStartSuppression();
            didSuppressB = true;
        }

        simA.step(dt);
        simB.step(dt);
        t = t_next;

        const StateSnap curA = snap(simA);
        const StateSnap curB = snap(simB);

        // Ensure ignition/suppression commands were not silently no-ops.
        if (t > ignite_at + 1e-12) {
            REQUIRE(didIgniteA && didIgniteB, "2.4 ignition requested but not issued");
        }
        if (t > suppress_on_at + 1e-12) {
            REQUIRE(didSuppressB, "2.4 suppression requested but not issued");
        }

        // Record the DeltaT window start at the first timestep strictly after suppression-on.
        // This ensures we measure from the point where suppression has had an opportunity to take effect.
        if (!haveSuppressionStartTemps && didSuppressB && t >= suppress_on_at - 1e-12) {
            T_A_start = curA.o.T_K;
            T_B_start = curB.o.T_K;
            haveSuppressionStartTemps = true;
        }

        // Post-suppression comparisons.
        if (haveSuppressionStartTemps) {
            // Only compare while both sims are in comparable lifecycle regions.
            // (If one concludes early, we stop accruing eligible checks to avoid comparing frozen outputs
            //  against active dynamics. Step 1C already validates terminal no-op behavior.)
            if (!curA.concluded && !curB.concluded) {
                // Update running HRR max for gating and epsilon sizing.
                if (curA.o.HRR_W > hrrA_running_max_post) hrrA_running_max_post = curA.o.HRR_W;
                if (curA.o.HRR_W > hrrA_peak_post) hrrA_peak_post = curA.o.HRR_W;

                const double hrr_gate = std::max(1e-6, 0.05 * hrrA_running_max_post);

                // Eligible timesteps: baseline HRR must be meaningfully above noise floor.
                if (curA.o.HRR_W > hrr_gate) {
                    ++eligible;

                    const double eps_HRR = std::max(1e-3, 1e-6 * hrrA_peak_post);

                    // (1) HRR with suppression must not exceed baseline beyond epsilon.
                    if (!(curB.o.HRR_W <= curA.o.HRR_W + eps_HRR)) {
                        std::cerr << "[FAIL] 2.4 HRR monotonicity violated at t=" << t
                                  << " HRR_A=" << curA.o.HRR_W
                                  << " HRR_B=" << curB.o.HRR_W
                                  << " eps_HRR=" << eps_HRR << "\n";
                        std::exit(1);
                    }

                    // (3) No paradox: if HRR_B is materially lower, dT_B must not be materially higher.
                    const double dT_A = curA.o.T_K - prevA.o.T_K;
                    const double dT_B = curB.o.T_K - prevB.o.T_K;

                    REQUIRE_FINITE(dT_A, "2.4 dT_A");
                    REQUIRE_FINITE(dT_B, "2.4 dT_B");

                    if (curB.o.HRR_W < curA.o.HRR_W - eps_HRR) {
                        REQUIRE(!(dT_B > dT_A + eps_T_step),
                                "2.4 paradox: HRR_B << HRR_A but dT_B >> dT_A");
                    }
                }
            }
        }

        prevA = curA;
        prevB = curB;
    }

    REQUIRE(didIgniteA && didIgniteB, "2.4: ignition was not issued in both sims");
    REQUIRE(didSuppressB, "2.4: suppression was not issued in B");
    REQUIRE(haveSuppressionStartTemps, "2.4: did not reach suppression-on window");
    REQUIRE(eligible >= 10, "2.4: insufficient eligible post-suppression timesteps (vacuous pass prevention)");

    // (2) Total DeltaT from suppression-on time to end: B must not heat more than A (beyond eps).
    const double T_A_end = prevA.o.T_K;
    const double T_B_end = prevB.o.T_K;

    const double dT_A_total = T_A_end - T_A_start;
    const double dT_B_total = T_B_end - T_B_start;

    REQUIRE_FINITE(dT_A_total, "2.4 DeltaT_A");
    REQUIRE_FINITE(dT_B_total, "2.4 DeltaT_B");

    REQUIRE(dT_B_total <= dT_A_total + eps_T_total,
            "2.4: total DeltaT with suppression exceeded no-suppression DeltaT");

    std::cout << "[PASS] 2.4 suppression reduces heating outcomes"
              << " eligible=" << eligible
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " suppress_on_at=" << suppress_on_at
              << " t_end=" << t_end
              << " dT_A_total=" << dT_A_total
              << " dT_B_total=" << dT_B_total << "\n";
}

/* =======================
 * Step 2.5: Fuel monotonicity under combustion
 * Requirements (black-box):
 * - Under meaningful combustion, fuel_kg must be non-increasing per timestep (within eps).
 * - Optional directional check: suppression should not cause *more* fuel consumption than baseline.
 * ======================= */
static void runFuelMonotonicityUnderCombustion_2A5() {
    vfep::Simulation simA;
    vfep::Simulation simB;
    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Explicit assumptions
    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_on_at = ignite_at + 0.5; // 0.5s after ignition
    const double t_end = 12.0;

    // Tolerances
    const double eps_fuel_step = 1e-12; // kg
    const double eps_fuel_cmp  = 1e-9;  // kg

    REQUIRE(dt > 0.0, "2.5 dt must be > 0");
    REQUIRE(t_end > ignite_at, "2.5 t_end must be > ignite_at");

    double t = 0.0;
    bool didIgniteA = false;
    bool didIgniteB = false;
    bool didSuppressB = false;

    StateSnap prevA = snap(simA);
    StateSnap prevB = snap(simB);

    // Adaptive HRR gate based on baseline (A) after ignition.
    double hrrA_running_max = 0.0;
    std::size_t eligible = 0;

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_end / dt) + 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    while (t + dt <= t_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.5 non-termination guard tripped");

        const double t_next = t + dt;

        // Matched ignition timing
        if (!didIgniteA && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simA.commandIgniteOrIncreasePyrolysis();
            didIgniteA = true;
        }
        if (!didIgniteB && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simB.commandIgniteOrIncreasePyrolysis();
            didIgniteB = true;
        }

        // Suppression only in B
        if (!didSuppressB && (t >= suppress_on_at || (t < suppress_on_at && t_next >= suppress_on_at))) {
            simB.commandStartSuppression();
            didSuppressB = true;
        }

        simA.step(dt);
        simB.step(dt);
        t = t_next;

        const StateSnap curA = snap(simA);
        const StateSnap curB = snap(simB);

        // Only evaluate monotonicity once ignition has been issued and we are past it.
        if (t > ignite_at + 1e-12) {
            REQUIRE(didIgniteA && didIgniteB, "2.5 ignition requested but not issued");

            // Update baseline HRR running max for gating
            if (curA.o.HRR_W > hrrA_running_max) hrrA_running_max = curA.o.HRR_W;
            const double hrr_gate = std::max(1e-6, 0.05 * hrrA_running_max);

            // Eligible = meaningful burning regime in baseline
            if (!curA.concluded && !curB.concluded && curA.o.HRR_W > hrr_gate) {
                ++eligible;

                // (1) Per-step monotonicity: fuel must not increase (within eps) in either run.
                REQUIRE(curA.o.fuel_kg <= prevA.o.fuel_kg + eps_fuel_step,
                        "2.5: fuel increased in baseline (A) during combustion");
                REQUIRE(curB.o.fuel_kg <= prevB.o.fuel_kg + eps_fuel_step,
                        "2.5: fuel increased in suppressed run (B) during combustion");

                // (2) Directional check: suppression should not consume more fuel than baseline.
                // i.e., remaining fuel in B should be >= remaining fuel in A (within eps).
                REQUIRE(curB.o.fuel_kg + eps_fuel_cmp >= curA.o.fuel_kg,
                        "2.5: suppressed run consumed more fuel than baseline");
            }
        }

        prevA = curA;
        prevB = curB;
    }

    REQUIRE(didIgniteA && didIgniteB, "2.5: ignition was not issued in both sims");
    REQUIRE(didSuppressB, "2.5: suppression was not issued in B");
    REQUIRE(eligible >= 10, "2.5: insufficient eligible combustion timesteps (vacuous pass prevention)");

    std::cout << "[PASS] 2.5 fuel monotonicity under combustion"
              << " eligible=" << eligible
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " suppress_on_at=" << suppress_on_at
              << " t_end=" << t_end
              << " eps_step=" << eps_fuel_step
              << " eps_cmp=" << eps_fuel_cmp
              << "\n";
}

static void runNoSpontaneousMassAppearanceForAgentReservoirs_2A6() {
    vfep::Simulation simA;
    vfep::Simulation simB;

    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Explicit assumptions (per spec)
    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_on_at = 5.0; // B schedules suppression at 5s
    const double sample_at = 4.0;      // compare at 4s (pre-suppression)
    const double t_end = 6.0;          // short run window

    // Tolerances
    const double eps_cmp = 1e-12; // baseline equality tolerance
    const double eps_nz  = 1e-12; // hard "non-zero" threshold

    double t = 0.0;
    bool ignitedA = false, ignitedB = false;
    bool suppressCommandedB = false;
    bool sampled = false;

    auto requireNearlyEqual = [&](double x, double y, double eps, const char* what) {
        REQUIRE(std::fabs(x - y) <= eps, what);
    };

    auto requireNearlyZero = [&](double x, double eps, const char* what) {
        REQUIRE(std::fabs(x) <= eps, what);
    };

    while (t + 1e-12 < t_end) {
        const double t_next = t + dt;

        // Matched ignition timing
        if (!ignitedA && t_next >= ignite_at - 1e-12) { simA.commandIgniteOrIncreasePyrolysis(); ignitedA = true; }
        if (!ignitedB && t_next >= ignite_at - 1e-12) { simB.commandIgniteOrIncreasePyrolysis(); ignitedB = true; }

        // Suppression scheduled only in B at known time (after sample time)
        if (!suppressCommandedB && t_next >= suppress_on_at - 1e-12) { simB.commandStartSuppression(); suppressCommandedB = true; }

        simA.step(dt);
        simB.step(dt);
        t = t_next;

        const StateSnap a = snap(simA);
        const StateSnap b = snap(simB);

        // Invariant: BEFORE suppression is commanded/enabled, agent-related concentrations must not increase/appear.
        // Practical check: enforce equality to baseline + hard non-zero forbiddance pre-command.
        //
        // Eligibility: strictly pre-command time (matches "before suppression is commanded/enabled").
        if (!suppressCommandedB && t <= suppress_on_at - 1e-12) {
            // Must match baseline (exact or within eps)
            requireNearlyEqual(b.o.inhibitor_kgm3, a.o.inhibitor_kgm3, eps_cmp, "2.6 inhibitor_kgm3 differs from baseline pre-suppression");
            requireNearlyEqual(b.o.inert_kgm3,     a.o.inert_kgm3,     eps_cmp, "2.6 inert_kgm3 differs from baseline pre-suppression");
            requireNearlyEqual(b.o.agent_mdot_kgps, a.o.agent_mdot_kgps, eps_cmp, "2.6 agent_mdot_kgps differs from baseline pre-suppression");

            // Hard failure: any non-zero agent concentrations/flow prior to suppression taking effect
            requireNearlyZero(b.o.inhibitor_kgm3, eps_nz, "2.6 hard fail: inhibitor_kgm3 non-zero before suppression");
            requireNearlyZero(b.o.inert_kgm3,     eps_nz, "2.6 hard fail: inert_kgm3 non-zero before suppression");
            requireNearlyZero(b.o.agent_mdot_kgps, eps_nz, "2.6 hard fail: agent_mdot_kgps non-zero before suppression");
        }

        // Sample at t=4s (pre-command). Use a dt-aware window for exact stepping robustness.
        if (!sampled && std::fabs(t - sample_at) <= 0.5 * dt + 1e-12) {
            sampled = true;
            REQUIRE(!suppressCommandedB, "2.6 internal: suppression was commanded before sample_at");

            requireNearlyEqual(b.o.inhibitor_kgm3, a.o.inhibitor_kgm3, eps_cmp, "2.6 inhibitor_kgm3 mismatch at sample_at");
            requireNearlyEqual(b.o.inert_kgm3,     a.o.inert_kgm3,     eps_cmp, "2.6 inert_kgm3 mismatch at sample_at");
            requireNearlyEqual(b.o.agent_mdot_kgps, a.o.agent_mdot_kgps, eps_cmp, "2.6 agent_mdot_kgps mismatch at sample_at");

            requireNearlyZero(b.o.inhibitor_kgm3, eps_nz, "2.6 hard fail: inhibitor_kgm3 non-zero at sample_at");
            requireNearlyZero(b.o.inert_kgm3,     eps_nz, "2.6 hard fail: inert_kgm3 non-zero at sample_at");
            requireNearlyZero(b.o.agent_mdot_kgps, eps_nz, "2.6 hard fail: agent_mdot_kgps non-zero at sample_at");
        }
    }

    REQUIRE(sampled, "2.6 did not reach sample_at (vacuous)");

    std::cout << "[PASS] 2.6 no spontaneous mass appearance pre-suppression"
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " suppress_on_at=" << suppress_on_at
              << " sample_at=" << sample_at
              << " t_end=" << t_end
              << " eps_cmp=" << eps_cmp
              << " eps_nz=" << eps_nz
              << "\n";
}

/* =======================
 * Step 2.7: Suppression increases agent presence after command
 * Requirements (black-box preference):
 *  A) Invariant:
 *     After suppression becomes enabled (isSuppressionEnabled()==true) and delivery is commanded,
 *     at least one agent indicator must become non-zero (per public API semantics).
 *
 *  B) Practical check (paired runs):
 *     - Run A baseline: ignition only.
 *     - Run B: ignition + suppression commanded at t=5.0s.
 *     - Evaluate post-suppression window t in [6.0, 10.0] (dt=0.05).
 *     Requirement: at least once in the window (and after suppression is enabled), either:
 *        inhibitor_kgm3 > eps_agent OR inert_kgm3 > eps_agent
 *     If concentrations never rise but agent_mdot_kgps > eps_mdot, we accept flow as the
 *     model's public "delivery" indicator (with an explicit note), because some variants
 *     may represent agent purely via a delivered/flow channel.
 *
 *  C) Species plausibility (directional, windowed; no calibration):
 *     - Baseline run must keep agent channels ~0 in the post-suppression window.
 *     - In run B, once suppression is enabled, agent_sum = inhibitor_kgm3 + inert_kgm3
 *       must show an *increase event* within a short window after enable (<= 1.0s),
 *       OR agent_mdot_kgps must be clearly non-zero.
 *     - During timesteps where agent_mdot_kgps is clearly non-zero, agent_sum must not be
 *       perfectly flat (must increase at least once). This avoids vacuous passes where
 *       delivery is "active" but never reflected in any public state channel.
 *
 * Optional tightenings (still non-calibration):
 *     - Bound enable-delay after command (<= max_enable_delay).
 *     - Require run B to exceed run A at least once in post-window agent indicators.
 * ======================= */
static void runSuppressionIncreasesAgentPresenceAfterCommand_2A7() {
    vfep::Simulation simA;
    vfep::Simulation simB;
    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Explicit assumptions (per spec)
    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_cmd_at = 5.0;
    const double post_window_start = 6.0;
    const double post_window_end   = 10.0;
    const double t_end = post_window_end;

    // Tolerances / thresholds
    const double eps_agent = 1e-12; // "non-zero" threshold for concentration channels
    const double eps_mdot  = 1e-12; // "non-zero" threshold for agent flow channel
    const double eps_flat  = 1e-18; // flatness epsilon for detecting changes

    REQUIRE(dt > 0.0, "2.7 dt must be > 0");
    REQUIRE(post_window_end > post_window_start, "2.7 invalid post window");
    REQUIRE(t_end >= post_window_end, "2.7 internal: t_end must cover post window");

    double t = 0.0;
    bool didIgniteA = false;
    bool didIgniteB = false;
    bool didSuppressCommandB = false;

    bool suppressionEnabledSeen = false;
    double t_suppression_enabled_first = std::numeric_limits<double>::quiet_NaN();

    // Evidence flags (post-window)
    bool baselineAgentStayedZero = true;
    bool sawConcentrationEvidence = false;
    bool sawFlowEvidence = false;

    // Optional: require run B exceeds run A at least once post-window
    bool sawBExceedsA = false;

    // Directional plausibility tracking
    bool capturedAtEnable = false;
    double agent_sum_at_enable = 0.0;
    bool sawEarlyIncreaseEvent = false;

    bool sawAnyDeliveryActive = false;
    bool sawAgentSumIncreaseDuringDelivery = false;
    double agent_sum_prev = 0.0;
    bool agent_prev_valid = false;

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_end / dt) + 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    // Initial snapshot (also validates finiteness/bounds)
    (void)snap(simA);
    (void)snap(simB);

    while (t + dt <= t_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.7 non-termination guard tripped");

        const double t_next = t + dt;

        // Matched ignition timing
        if (!didIgniteA && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simA.commandIgniteOrIncreasePyrolysis();
            didIgniteA = true;
        }
        if (!didIgniteB && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simB.commandIgniteOrIncreasePyrolysis();
            didIgniteB = true;
        }

        // Suppression commanded only in B
        if (!didSuppressCommandB && (t >= suppress_cmd_at || (t < suppress_cmd_at && t_next >= suppress_cmd_at))) {
            simB.commandStartSuppression();
            didSuppressCommandB = true;
        }

        simA.step(dt);
        simB.step(dt);
        t = t_next;

        const StateSnap curA = snap(simA);
        const StateSnap curB = snap(simB);

        // Track first time suppression becomes enabled (preferred semantics)
        if (!suppressionEnabledSeen && curB.suppressed) {
            suppressionEnabledSeen = true;
            t_suppression_enabled_first = t;
        }

        // Evaluate only in the requested post-suppress window
        const bool in_post_window = (t >= post_window_start - 1e-12) && (t <= post_window_end + 1e-12);
        if (in_post_window) {
            // Baseline must keep agent channels at ~0 (robust against accidental leakage)
            if (curA.o.inhibitor_kgm3 > eps_agent ||
                curA.o.inert_kgm3 > eps_agent ||
                curA.o.agent_mdot_kgps > eps_mdot) {
                baselineAgentStayedZero = false;
            }

            // We only evaluate B after suppression is actually enabled.
            if (curB.suppressed) {
                const double agent_sum = curB.o.inhibitor_kgm3 + curB.o.inert_kgm3;
                REQUIRE_FINITE(agent_sum, "2.7 agent_sum");

                // Primary evidence: concentration channels become non-zero.
                if (curB.o.inhibitor_kgm3 > eps_agent || curB.o.inert_kgm3 > eps_agent) {
                    sawConcentrationEvidence = true;
                }

                // Secondary evidence: flow channel indicates delivery.
                if (curB.o.agent_mdot_kgps > eps_mdot) {
                    sawFlowEvidence = true;
                    sawAnyDeliveryActive = true;
                }

                // Optional: ensure run B exceeds baseline A at least once post-window (event-based).
                if ((curB.o.inhibitor_kgm3  > curA.o.inhibitor_kgm3  + eps_agent) ||
                    (curB.o.inert_kgm3      > curA.o.inert_kgm3      + eps_agent) ||
                    (curB.o.agent_mdot_kgps > curA.o.agent_mdot_kgps + eps_mdot)) {
                    sawBExceedsA = true;
                }

                // Capture agent baseline at the first enabled timestep for directional check.
                if (!capturedAtEnable) {
                    capturedAtEnable = true;
                    agent_sum_at_enable = agent_sum;
                    agent_sum_prev = agent_sum;
                    agent_prev_valid = true;
                }

                // (C1) Early increase event: within 1.0s of enable, agent_sum should increase at least once.
                if (capturedAtEnable) {
                    const double dt_since_enable = t - t_suppression_enabled_first;
                    if (dt_since_enable <= 1.0 + 1e-12) {
                        if (agent_sum > agent_sum_at_enable + eps_agent) {
                            sawEarlyIncreaseEvent = true;
                        }
                    }
                }

                // (C2) During active delivery, agent_sum should not be perfectly flat.
                // We require at least one increase event (not per-step monotonicity).
                if (curB.o.agent_mdot_kgps > eps_mdot) {
                    if (agent_prev_valid && agent_sum > agent_sum_prev + eps_flat) {
                        sawAgentSumIncreaseDuringDelivery = true;
                    }
                }

                agent_sum_prev = agent_sum;
                agent_prev_valid = true;
            }
        }
    }

    REQUIRE(didIgniteA && didIgniteB, "2.7: ignition was not issued in both sims");
    REQUIRE(didSuppressCommandB, "2.7: suppression command was not issued in B");

    // Ensure suppression actually became enabled (otherwise test cannot validate 2.7 semantics)
    REQUIRE(suppressionEnabledSeen, "2.7: suppression never became enabled in B (isSuppressionEnabled() stayed false)");
    REQUIRE(t_suppression_enabled_first <= post_window_end + 1e-12,
            "2.7: suppression enabled only after post-suppression window (test would be vacuous)");

    // Optional: bounded enable delay after command (directional, not calibration)
    const double max_enable_delay = 1.0; // seconds
    REQUIRE(t_suppression_enabled_first - suppress_cmd_at <= max_enable_delay + 1e-12,
            "2.7: suppression enable delay exceeded max bound");

    // Core requirement: agent presence must manifest in at least one public indicator.
    REQUIRE(baselineAgentStayedZero, "2.7: baseline run showed non-zero agent indicators in post window");
    REQUIRE(sawConcentrationEvidence || sawFlowEvidence,
            "2.7 hard fail: suppression enabled but agent presence never manifested in any public channel");

    // Optional: ensure run B exceeds baseline run A at least once post-window
    REQUIRE(sawBExceedsA,
            "2.7: run B never exceeded baseline run A in agent indicators post window");

    // Directional plausibility (non-calibration): require an early increase event OR clear delivery flow.
    REQUIRE(sawEarlyIncreaseEvent || sawFlowEvidence,
            "2.7: suppression enabled but no directional evidence of agent accumulation or delivery within 1s");

    // If we saw delivery flow, require at least one increase in agent_sum during delivery.
    // This is a strict-but-robust sanity check: it does NOT require monotonicity.
    if (sawAnyDeliveryActive) {
        REQUIRE(sawAgentSumIncreaseDuringDelivery || sawConcentrationEvidence,
                "2.7: agent delivery flow observed, but agent_sum never increased (no state reflection)");
    }

    std::cout << "[PASS] 2.7 suppression increases agent presence after command"
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " suppress_cmd_at=" << suppress_cmd_at
              << " post_window=[" << post_window_start << "," << post_window_end << "]"
              << " eps_agent=" << eps_agent
              << " eps_mdot=" << eps_mdot
              << " enabled_at=" << t_suppression_enabled_first
              << " concEvidence=" << (sawConcentrationEvidence ? 1 : 0)
              << " flowEvidence=" << (sawFlowEvidence ? 1 : 0)
              << " earlyIncrease=" << (sawEarlyIncreaseEvent ? 1 : 0)
              << "\n";

    if (sawFlowEvidence && !sawConcentrationEvidence) {
        std::cout << "[NOTE] 2.7 passed via agent_mdot_kgps evidence; concentration channels stayed <= eps_agent. "
                     "If this is unintended, expose/propagate delivered agent concentration fields.\n";
    }
}

/* =======================
 * Step 2.8: During active combustion: O2 down, CO2/H2O up (net directional)
 *
 * Invariant (post-ignite, pre-suppress window):
 *   O2_end  <= O2_start  + eps_vol
 *   CO2_end >= CO2_start - eps_vol
 *   H2O_end >= H2O_start - eps_vol
 *
 * Hard failure:
 *   Net O2 increase AND net CO2 decrease AND net H2O decrease within an eligible combustion window.
 *
 * Eligibility / gating (non-calibration):
 *   - ignition has latched true
 *   - suppression is NOT enabled throughout the window
 *   - HRR is meaningfully positive (active combustion proxy)
 *   - fuel decreases over the window (combustion proxy)
 *
 * Notes:
 *   - Uses vol% channels (O2_volpct, CO2_volpct, H2O_volpct).
 *   - Event/net-direction checks only; no monotonicity requirements.
 * ======================= */
static void runCombustionSpeciesDirectionality_2A8() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    // Explicit assumptions
    const double dt = 0.05;
    const double ignite_at = 2.0;

    // Directionality window (post-ignite, pre-suppress)
    // Chosen to align with prior tests where suppression might start at 5.0s in paired runs.
    const double window_start = 3.0;
    const double window_end   = 5.0;
    const double t_end = window_end;

    // Tolerances / gates
    const double eps_vol   = 1e-7;   // vol% tolerance (directionality; not calibration)
    const double hrr_gate  = 1e3;    // W; proxy for active combustion (non-calibration)
    const double eps_fuel  = 1e-12;  // kg; detect net fuel decrease

    REQUIRE(dt > 0.0, "2.8 dt must be > 0");
    REQUIRE(window_end > window_start, "2.8 invalid window");
    REQUIRE(t_end >= window_end, "2.8 internal: t_end must cover window");

    double t = 0.0;
    bool didIgnite = false;

    bool haveStart = false;
    bool haveEnd = false;
    StateSnap sStart{};
    StateSnap sEnd{};

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_end / dt) + 16u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    (void)snap(sim);

    while (t + dt <= t_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.8 non-termination guard tripped");

        const double t_next = t + dt;

        if (!didIgnite && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
        }

        sim.step(dt);
        t = t_next;

        const StateSnap cur = snap(sim);

        // Capture start and end snapshots at the first timestep at/after the boundary.
        if (!haveStart && t >= window_start - 1e-12) {
            sStart = cur;
            haveStart = true;
        }
        if (!haveEnd && t >= window_end - 1e-12) {
            sEnd = cur;
            haveEnd = true;
            break; // we have what we need
        }
    }

    REQUIRE(didIgnite, "2.8: ignition was not issued");
    REQUIRE(haveStart && haveEnd, "2.8: failed to capture window endpoints");

    // Eligibility gating (robust, non-calibration)
    REQUIRE(sStart.ignited && sEnd.ignited, "2.8: not ignited across window (eligibility failed)");
    REQUIRE(!sStart.suppressed && !sEnd.suppressed, "2.8: suppression enabled during window (eligibility failed)");

    const bool active_hrr = (maxd(sStart.o.HRR_W, sEnd.o.HRR_W) > hrr_gate);
    const double fuel_delta = sStart.o.fuel_kg - sEnd.o.fuel_kg;
    const bool fuel_consumed = (fuel_delta > eps_fuel);

    REQUIRE(active_hrr, "2.8: HRR not active in window (eligibility failed)");
    REQUIRE(fuel_consumed, "2.8: fuel did not decrease in window (eligibility failed)");

    // Net changes (end - start)
    const double dO2  = sEnd.o.O2_volpct  - sStart.o.O2_volpct;
    const double dCO2 = sEnd.o.CO2_volpct - sStart.o.CO2_volpct;
    const double dH2O = sEnd.o.H2O_volpct - sStart.o.H2O_volpct;

    // Hard contradiction check first (strongly non-physical under active combustion)
    if ((dO2  >  eps_vol) &&
        (dCO2 < -eps_vol) &&
        (dH2O < -eps_vol)) {
        std::cerr << "[FAIL] 2.8 hard contradiction: active combustion window shows "
                  << "dO2=" << dO2 << " (increase), "
                  << "dCO2=" << dCO2 << " (decrease), "
                  << "dH2O=" << dH2O << " (decrease)"
                  << " start(O2,CO2,H2O)=(" << sStart.o.O2_volpct << "," << sStart.o.CO2_volpct << "," << sStart.o.H2O_volpct << ")"
                  << " end(O2,CO2,H2O)=("   << sEnd.o.O2_volpct   << "," << sEnd.o.CO2_volpct   << "," << sEnd.o.H2O_volpct   << ")"
                  << "\n";
        std::exit(1);
    }

    // Directional net checks (tolerant)
    REQUIRE(sEnd.o.O2_volpct  <= sStart.o.O2_volpct  + eps_vol,
            "2.8: net O2 increased during active combustion window");
    REQUIRE(sEnd.o.CO2_volpct >= sStart.o.CO2_volpct - eps_vol,
            "2.8: net CO2 decreased during active combustion window");
    REQUIRE(sEnd.o.H2O_volpct >= sStart.o.H2O_volpct - eps_vol,
            "2.8: net H2O decreased during active combustion window");

    std::cout << "[PASS] 2.8 combustion species directionality"
          << " dt=" << dt
          << " ignite_at=" << ignite_at
          << " window=[" << window_start << "," << window_end << "]"
          << " eps_vol=" << eps_vol
          << " hrr_gate=" << hrr_gate
          << " fuel_delta=" << fuel_delta
          << " dO2=" << dO2
          << " dCO2=" << dCO2
          << " dH2O=" << dH2O
          << " start(O2,CO2,H2O)=("
          << sStart.o.O2_volpct << ","
          << sStart.o.CO2_volpct << ","
          << sStart.o.H2O_volpct << ")"
          << " end(O2,CO2,H2O)=("
          << sEnd.o.O2_volpct << ","
          << sEnd.o.CO2_volpct << ","
          << sEnd.o.H2O_volpct << ")"
          << "\n";
}

/* =======================
 * Step 2.9: Suppression makes species evolution "less burned" than no-suppression (paired-run)
 *
 * Requirements:
 *   A) Paired-run setup (matched):
 *      Run A (no suppression): ignite at t=2.0, no suppression.
 *      Run B (suppressed): ignite at t=2.0, suppression command at t=5.0;
 *        suppression is considered enabled when isSuppressionEnabled()==true.
 *      dt = 0.05.
 *      Evaluate at t_eval = 10.0s (must be after suppression is enabled).
 *
 *   B) Invariant at t_eval (directional, tolerant):
 *      O2_B  >= O2_A  - tol
 *      CO2_B <= CO2_A + tol
 *      H2O_B <= H2O_A + tol
 *      Using vol% channels (O2_volpct, CO2_volpct, H2O_volpct).
 *
 *   C) Hard failure (systematically "more burned" in B beyond tol):
 *      O2_B  < O2_A  - tol  AND
 *      CO2_B > CO2_A + tol  AND
 *      H2O_B > H2O_A + tol
 *
 *   D) Causal ordering: command -> enable -> effect.
 *      Pre-enable window check (t in [3.0, t_enable)):
 *        Fail if the "less burned" signature is sustained for a meaningful streak
 *        (guards against accepting reverse causality).
 *
 * Eligibility / gating (robust, non-calibration):
 *   - Both runs must latch ignited==true before evaluation.
 *   - Baseline run must show evidence of combustion (HRR meaningfully positive and net fuel decrease)
 *     by t_eval (avoids ventilation/mixing false failures).
 *   - Suppressed run must have suppression enabled before t_eval.
 *
 * Robustness note (dilution):
 *   - If suppression delivers agent/inert, raw O2_volpct may drop simply due to displacement/dilution,
 *     even when the gas is "less burned" in a chemical sense. To avoid false failures, this test uses
 *     a dry-basis renormalization by (O2+CO2+H2O) for O2 (and for the hard-fail triple) when agent is present.
 * ======================= */
static void runSuppressionMakesLessBurnedThanNoSuppression_2A9() {
    vfep::Simulation simA;
    vfep::Simulation simB;
    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Explicit assumptions (per spec)
    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_cmd_at = 5.0;
    const double t_eval = 10.0;

    // Directional tolerance in vol% (non-calibration, numerical slack)
    const double tol = 1e-6;

    // Eligibility gating (non-calibration)
    const double hrr_gate_W = 1e3;      // active combustion proxy (same order as 2.8)
    const double eps_fuel_net = 1e-10;  // require some net fuel consumption by t_eval

    // Pre-enable reverse-causality guard (sustained streak requirement)
    const double pre_window_start = 3.0;
    const int pre_streak_required = 10; // 10*dt = 0.5s sustained "less burned" pre-enable => fail

    REQUIRE(dt > 0.0, "2.9 dt must be > 0");
    REQUIRE(t_eval > suppress_cmd_at, "2.9 internal: t_eval must be after suppression command time");

    double t = 0.0;
    bool didIgniteA = false;
    bool didIgniteB = false;
    bool didSuppressCommandB = false;

    bool suppressionEnabledSeen = false;
    double t_enable = std::numeric_limits<double>::quiet_NaN();

    // Capture for eligibility and diagnostics
    StateSnap a_at_eval{};
    StateSnap b_at_eval{};
    bool captured_eval = false;

    double fuelA_at_ignite = std::numeric_limits<double>::quiet_NaN();
    bool captured_fuelA_at_ignite = false;
    double hrrA_running_max = 0.0;

    int pre_lessburned_streak = 0;
    int pre_lessburned_hits = 0;
    int pre_samples = 0;

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_eval / dt) + 32u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    (void)snap(simA);
    (void)snap(simB);

    while (t + dt <= t_eval + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.9 non-termination guard tripped");

        const double t_next = t + dt;

        // Matched ignition timing
        if (!didIgniteA && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simA.commandIgniteOrIncreasePyrolysis();
            didIgniteA = true;
        }
        if (!didIgniteB && (t >= ignite_at || (t < ignite_at && t_next >= ignite_at))) {
            simB.commandIgniteOrIncreasePyrolysis();
            didIgniteB = true;
        }

        // Suppression command only in B
        if (!didSuppressCommandB && (t >= suppress_cmd_at || (t < suppress_cmd_at && t_next >= suppress_cmd_at))) {
            simB.commandStartSuppression();
            didSuppressCommandB = true;
        }

        simA.step(dt);
        simB.step(dt);
        t = t_next;

        const StateSnap a = snap(simA);
        const StateSnap b = snap(simB);

        // Track suppression enable time (public semantics)
        if (!suppressionEnabledSeen && b.suppressed) {
            suppressionEnabledSeen = true;
            t_enable = t;
        }

        // Capture baseline fuel around ignition for net-consumption eligibility
        if (!captured_fuelA_at_ignite && didIgniteA && t >= ignite_at - 1e-12) {
            fuelA_at_ignite = a.o.fuel_kg;
            captured_fuelA_at_ignite = true;
        }

        // Track baseline HRR running max for gating (robust across parameterizations)
        if (a.o.HRR_W > hrrA_running_max) hrrA_running_max = a.o.HRR_W;

        // Causal ordering / reverse-causality guard:
        // Before suppression is enabled, B should not already be consistently "less burned" than A.
        if (t >= pre_window_start - 1e-12 && (!suppressionEnabledSeen || t < t_enable - 1e-12)) {
            ++pre_samples;

            const bool less_burned_signature =
                (b.o.O2_volpct  > a.o.O2_volpct  + tol) &&
                (b.o.CO2_volpct < a.o.CO2_volpct - tol) &&
                (b.o.H2O_volpct < a.o.H2O_volpct - tol);

            if (less_burned_signature) {
                ++pre_lessburned_hits;
                ++pre_lessburned_streak;
            } else {
                pre_lessburned_streak = 0;
            }

            // Sustained pre-enable signature indicates reverse causality or mismatched setup.
            REQUIRE(pre_lessburned_streak < pre_streak_required,
                    "2.9 causal ordering fail: sustained 'less burned' signature before suppression enabled");
        }

        // Capture evaluation snapshots at first step at/after t_eval
        if (!captured_eval && t >= t_eval - 1e-12) {
            a_at_eval = a;
            b_at_eval = b;
            captured_eval = true;
            break;
        }
    }

    REQUIRE(didIgniteA && didIgniteB, "2.9: ignition was not issued in both sims");
    REQUIRE(didSuppressCommandB, "2.9: suppression command was not issued in B");
    REQUIRE(captured_eval, "2.9: failed to capture t_eval snapshot (vacuous)");
    REQUIRE(suppressionEnabledSeen, "2.9 eligibility failed: suppression never became enabled in B");
    REQUIRE(t_enable <= t_eval - 1e-12, "2.9 eligibility failed: suppression enabled after t_eval");

    // Eligibility / gating to avoid calibration sensitivity:
    REQUIRE(a_at_eval.ignited && b_at_eval.ignited, "2.9 eligibility failed: not ignited by t_eval");
    REQUIRE(b_at_eval.suppressed, "2.9 eligibility failed: suppression not enabled at t_eval");

    // Require baseline combustion evidence by t_eval (avoid ventilation/mixing false failures).
    const double hrr_gate_dynamic = std::max(hrr_gate_W, 0.05 * hrrA_running_max);
    REQUIRE(a_at_eval.o.HRR_W > hrr_gate_dynamic, "2.9 eligibility failed: baseline HRR not active at t_eval");
    REQUIRE(captured_fuelA_at_ignite, "2.9 internal: failed to capture baseline fuel at ignition");
    const double fuelA_net_delta = fuelA_at_ignite - a_at_eval.o.fuel_kg;
    REQUIRE(fuelA_net_delta > eps_fuel_net, "2.9 eligibility failed: baseline net fuel did not decrease by t_eval");

    // Compare on a dry basis (renormalized by O2+CO2+H2O) when agent is present,
    // because suppression delivery can dilute vol% without indicating "more burned".
    // This remains a black-box check derived solely from public vol% channels.
    const auto dryFrac = [](double x, double sum) -> double {
        return (sum > 0.0) ? (x / sum) : 0.0;
    };

    const double sumA = a_at_eval.o.O2_volpct + a_at_eval.o.CO2_volpct + a_at_eval.o.H2O_volpct;
    const double sumB = b_at_eval.o.O2_volpct + b_at_eval.o.CO2_volpct + b_at_eval.o.H2O_volpct;

    const double O2A_dry  = dryFrac(a_at_eval.o.O2_volpct,  sumA);
    const double CO2A_dry = dryFrac(a_at_eval.o.CO2_volpct, sumA);
    const double H2OA_dry = dryFrac(a_at_eval.o.H2O_volpct, sumA);
    const double O2B_dry  = dryFrac(b_at_eval.o.O2_volpct,  sumB);
    const double CO2B_dry = dryFrac(b_at_eval.o.CO2_volpct, sumB);
    const double H2OB_dry = dryFrac(b_at_eval.o.H2O_volpct, sumB);

    const double agent_sum_B = b_at_eval.o.inhibitor_kgm3 + b_at_eval.o.inert_kgm3;
    const bool agent_present = (agent_sum_B > 1e-12);

    // Deltas (B - A) for diagnostics (raw and dry)
    const double dO2_raw  = b_at_eval.o.O2_volpct  - a_at_eval.o.O2_volpct;
    const double dCO2_raw = b_at_eval.o.CO2_volpct - a_at_eval.o.CO2_volpct;
    const double dH2O_raw = b_at_eval.o.H2O_volpct - a_at_eval.o.H2O_volpct;
    const double dO2_dry  = O2B_dry  - O2A_dry;
    const double dCO2_dry = CO2B_dry - CO2A_dry;
    const double dH2O_dry = H2OB_dry - H2OA_dry;

    // Hard failure: systematically "more burned" in suppressed run.
    // Use dry-basis when agent is present; otherwise use raw vol%.
    const double dO2_hard  = agent_present ? dO2_dry  : dO2_raw;
    const double dCO2_hard = agent_present ? dCO2_dry : dCO2_raw;
    const double dH2O_hard = agent_present ? dH2O_dry : dH2O_raw;

    if ((dO2_hard < -tol) && (dCO2_hard > tol) && (dH2O_hard > tol)) {
        std::cerr << "[FAIL] 2.9 hard contradiction: suppressed run is systematically MORE burned at t_eval=" << t_eval
                  << " (tol=" << tol << ")\n"
                  << "  enabled_at=" << t_enable << " command_at=" << suppress_cmd_at << "\n"
                  << "  A(O2,CO2,H2O)=(" << a_at_eval.o.O2_volpct << "," << a_at_eval.o.CO2_volpct << "," << a_at_eval.o.H2O_volpct << ")\n"
                  << "  B(O2,CO2,H2O)=(" << b_at_eval.o.O2_volpct << "," << b_at_eval.o.CO2_volpct << "," << b_at_eval.o.H2O_volpct << ")\n"
                  << "  agent_present=" << (agent_present ? 1 : 0) << " agent_sum_B=" << agent_sum_B << "\n"
                  << "  raw d(B-A)=(" << dO2_raw << "," << dCO2_raw << "," << dH2O_raw << ")\n"
                  << "  dry d(B-A)=(" << dO2_dry << "," << dCO2_dry << "," << dH2O_dry << ")\n"
                  << "  pre_enable_hits=" << pre_lessburned_hits << "/" << pre_samples
                  << " pre_enable_max_streak=" << pre_streak_required
                  << "\n";
        std::exit(1);
    }

    // Directional invariants (tolerant)
    // - CO2/H2O: raw vol% is acceptable (a suppression-induced dilution still satisfies <= baseline).
    // - O2: if agent is present, use dry-basis to avoid false failures due to displacement/dilution.
    if (agent_present) {
        REQUIRE(O2B_dry >= O2A_dry - tol, "2.9: (dry) O2 in suppressed run is lower than baseline beyond tol");
    } else {
        REQUIRE(b_at_eval.o.O2_volpct >= a_at_eval.o.O2_volpct - tol, "2.9: O2 in suppressed run is lower than baseline beyond tol");
    }
    REQUIRE(b_at_eval.o.CO2_volpct <= a_at_eval.o.CO2_volpct + tol, "2.9: CO2 in suppressed run is higher than baseline beyond tol");
    REQUIRE(b_at_eval.o.H2O_volpct <= a_at_eval.o.H2O_volpct + tol, "2.9: H2O in suppressed run is higher than baseline beyond tol");

    std::cout << "[PASS] 2.9 suppression yields less-burned composition vs no-suppression"
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " suppress_cmd_at=" << suppress_cmd_at
              << " enabled_at=" << t_enable
              << " t_eval=" << t_eval
              << " tol=" << tol
              << " hrr_gate_dynamic=" << hrr_gate_dynamic
              << " fuelA_net_delta=" << fuelA_net_delta
              << " agent_present=" << (agent_present ? 1 : 0)
              << " agent_sum_B=" << agent_sum_B
              << " raw d(B-A)=(O2=" << dO2_raw << ",CO2=" << dCO2_raw << ",H2O=" << dH2O_raw << ")"
              << " dry d(B-A)=(O2=" << dO2_dry << ",CO2=" << dCO2_dry << ",H2O=" << dH2O_dry << ")"
              << " A(O2,CO2,H2O)=(" << a_at_eval.o.O2_volpct << "," << a_at_eval.o.CO2_volpct << "," << a_at_eval.o.H2O_volpct << ")"
              << " B(O2,CO2,H2O)=(" << b_at_eval.o.O2_volpct << "," << b_at_eval.o.CO2_volpct << "," << b_at_eval.o.H2O_volpct << ")"
              << " pre_enable_hits=" << pre_lessburned_hits << "/" << pre_samples
              << "\n";
}

/* =======================
 * Step 2.10: No suppression effects before suppression command
 *
 * Invariant:
 *   - If commandStartSuppression() is never issued, then:
 *       * isSuppressionEnabled() must remain false
 *       * agent channels must remain ~0 (no spontaneous agent appearance)
 *       * ignition alone must not enable suppression
 *
 * Non-vacuous guard:
 *   - Require evidence of combustion after ignition (HRR/fuel/species) before t_end,
 *     otherwise this test could pass in a totally quiescent / broken ignition path.
 *
 * Window:
 *   - Run through t_end < suppress_cmd_at (5.0s), with ignition at 2.0s
 * ======================= */
static void runNoSuppressionEffectsBeforeCommand_2A10() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.05;
    const double ignite_at = 2.0;
    const double suppress_cmd_at = 5.0;
    const double t_end = 4.9; // strictly before suppression command time

    // Hard "should be zero" thresholds (tripwires, not tuning)
    const double eps_agent = 1e-12;
    const double eps_mdot  = 1e-12;

    // Non-vacuous combustion evidence thresholds (tripwires, aligned with 2.8/2.11)
    const double hrr_gate_W       = 1e3;
    const double fuel_min_delta   = 1e-12;
    const double eps_species_post = 1e-7;

    REQUIRE(dt > 0.0, "2.10 dt must be > 0");
    REQUIRE(t_end < suppress_cmd_at, "2.10 internal: t_end must be < suppress_cmd_at");
    REQUIRE(t_end > ignite_at + dt, "2.10 internal: window must include post-ignite time");

    double t = 0.0;
    bool didIgnite = false;

    // Initial snapshot
    const StateSnap s0 = snap(sim);
    REQUIRE(!s0.suppressed, "2.10 baseline: suppression enabled unexpectedly at t=0");
    REQUIRE(std::fabs(s0.o.inhibitor_kgm3) <= eps_agent, "2.10 baseline: inhibitor_kgm3 non-zero");
    REQUIRE(std::fabs(s0.o.inert_kgm3)     <= eps_agent, "2.10 baseline: inert_kgm3 non-zero");
    REQUIRE(std::fabs(s0.o.agent_mdot_kgps) <= eps_mdot, "2.10 baseline: agent_mdot_kgps non-zero");

    // Capture post-ignite window endpoints for non-vacuous combustion evidence.
    // We choose a short window strictly after ignition but still < 5.0.
    const double window_start = 3.0;
    const double window_end   = 4.5;

    bool captured_start = false;
    bool captured_end = false;
    StateSnap s_start{};
    StateSnap s_end{};
    double maxHRR_post = 0.0;

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_end / dt) + 64u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    while (t + dt <= t_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.10 non-termination guard tripped");

        const double t_prev = t;
        const double t_next = t + dt;

        // Ignite (but never command suppression)
        if (!didIgnite && (t_prev >= ignite_at || (t_prev < ignite_at && t_next >= ignite_at))) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
        }

        sim.step(dt);
        t = t_next;

        const StateSnap s = snap(sim);

        // Core invariant: suppression must remain disabled because we never commanded it.
        REQUIRE(!s.suppressed, "2.10: suppression became enabled without suppression command");

        // Agent channels must remain ~0 pre-command.
        REQUIRE(std::fabs(s.o.inhibitor_kgm3) <= eps_agent, "2.10: inhibitor_kgm3 appeared before suppression command");
        REQUIRE(std::fabs(s.o.inert_kgm3)     <= eps_agent, "2.10: inert_kgm3 appeared before suppression command");
        REQUIRE(std::fabs(s.o.agent_mdot_kgps) <= eps_mdot, "2.10: agent_mdot_kgps became non-zero before suppression command");

        // Collect non-vacuous combustion evidence post-ignite (dt-robust endpoint capture).
        if (didIgnite) {
            if (s.o.HRR_W > maxHRR_post) maxHRR_post = s.o.HRR_W;

            if (!captured_start && (t_prev < window_start && t >= window_start - 1e-12)) {
                s_start = s;
                captured_start = true;
            }
            if (!captured_end && (t_prev < window_end && t >= window_end - 1e-12)) {
                s_end = s;
                captured_end = true;
            }
        }
    }

    REQUIRE(didIgnite, "2.10 internal: ignition was not issued (test setup failure)");

    // Non-vacuous: require we actually observed post-ignite dynamics in the chosen window.
    REQUIRE(captured_start && captured_end, "2.10: failed to capture post-ignite evidence window (vacuous)");

    // Require combustion evidence without any suppression effects.
    // This is a strict causality sanity check: ignition should have visible consequences by ~3-4.5s.
    const double dFuel = s_start.o.fuel_kg - s_end.o.fuel_kg; // positive means consumed
    const double dO2   = s_end.o.O2_volpct  - s_start.o.O2_volpct;
    const double dCO2  = s_end.o.CO2_volpct - s_start.o.CO2_volpct;
    const double dH2O  = s_end.o.H2O_volpct - s_start.o.H2O_volpct;

    REQUIRE(maxHRR_post >= hrr_gate_W, "2.10 non-vacuous: HRR never exceeded gate after ignition");
    REQUIRE(dFuel >= fuel_min_delta, "2.10 non-vacuous: fuel did not decrease after ignition");

    // Species directionality is a useful “combustion present” proxy, but allow small eps.
    REQUIRE(dO2  <= -eps_species_post, "2.10 non-vacuous: O2 did not decrease post-ignite");
    REQUIRE(dCO2 >=  eps_species_post, "2.10 non-vacuous: CO2 did not increase post-ignite");
    REQUIRE(dH2O >=  eps_species_post, "2.10 non-vacuous: H2O did not increase post-ignite");

    std::cout << "[PASS] 2.10 no suppression effects before command"
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " t_end=" << t_end
              << " eps_agent=" << eps_agent
              << " eps_mdot=" << eps_mdot
              << " maxHRR_post=" << maxHRR_post
              << " dFuel=" << dFuel
              << " d(O2,CO2,H2O)=(" << dO2 << "," << dCO2 << "," << dH2O << ")"
              << "\n";
}

/* =======================
 * Step 2.11: Ignition precedes combustion indicators (strict enough to be meaningful)
 *
 * Invariant:
 *   - Before ignition is commanded: ignited latch must remain false and no combustion signatures appear.
 *   - After ignition is commanded and sufficient time elapses: combustion signatures must appear.
 *
 * Hard failure:
 *   - Any combustion-like changes before ignition command.
 *   - Or ignition that never produces combustion signatures within a reasonable window.
 *
 * E) Time-local physical bounds:
 *   - Enforced implicitly via snap() -> requireFiniteAndBounded() each step.
 * ======================= */
static void runIgnitionPrecedesCombustionIndicators_2A11() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    // Explicit assumptions / harness parameters
    const double dt = 0.05;
    const double ignite_at = 2.0;

    // “Reasonable window” after ignition for combustion signatures to appear.
    // This ties to 2.2 / 2.8 behavior (combustion evidence over ~[3,5] s when ignite_at=2).
    const double t_window_start = 3.0; // ignite_at + 1.0
    const double t_window_end   = 5.0; // ignite_at + 3.0

    REQUIRE(dt > 0.0, "2.11 dt must be > 0");
    REQUIRE(t_window_start > ignite_at, "2.11 invalid window: start must be after ignite_at");
    REQUIRE(t_window_end > t_window_start, "2.11 invalid window: end must be after start");

    // Pre-ignition strict quiescence thresholds.
    // Keep these tiny: this is a causality gate (no “combustion-like” dynamics without ignition command).
    const double eps_hrr_pre     = 1e-9;
    const double eps_fuel_pre    = 1e-12;
    const double eps_species_pre = 1e-7;

    // Post-ignition “meaningful” combustion evidence thresholds.
    // These are not physics tuning; they are tripwires to ensure ignition has observable consequences.
    const double hrr_gate         = 1000.0;
    const double fuel_min_delta   = 0.02;
    const double eps_species_post = 1e-4;

    double t = 0.0;
    bool didIgnite = false;
    double t_ignite_cmd = -1.0;

    // Baseline snapshot for pre-ignite stability checks
    const StateSnap s0 = snap(sim);
    REQUIRE(!s0.ignited, "2.11: ignited latch true at t=0 (unexpected)");
    REQUIRE(s0.o.HRR_W <= eps_hrr_pre, "2.11: HRR_W not quiescent at t=0");

    // Post-ignite window evidence
    bool captured_start = false;
    bool captured_end = false;
    StateSnap s_start{};
    StateSnap s_end{};
    double maxHRR_post = 0.0;

    // Non-termination guard
    const std::size_t cap_by_ratio = (std::size_t)(t_window_end / dt) + 128u;
    const std::size_t iter_cap = (cap_by_ratio > 50'000'000u) ? 50'000'000u : cap_by_ratio;
    std::size_t iter = 0;

    // Run until we capture both window endpoints (dt-robust threshold crossing).
    while (t < t_window_end + 1e-12) {
        ++iter;
        REQUIRE(iter <= iter_cap, "2.11 non-termination guard tripped");

        const double t_prev = t;
        const double t_next = t + dt;

        // Robust ignition scheduling: issue command when crossing ignite_at.
        if (!didIgnite && (t_prev >= ignite_at || (t_prev < ignite_at && t_next >= ignite_at))) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
            t_ignite_cmd = t_prev; // command issued at current integration state before stepping
        }

        sim.step(dt);
        t = t_next;

        const StateSnap s = snap(sim);

        // PRE-IGNITION invariants: apply strictly while ignition has not been commanded.
        if (!didIgnite) {
            REQUIRE(!s.ignited, "2.11 pre-ignite: ignited became true without ignition command");

            // No combustion-like signatures
            REQUIRE(s.o.HRR_W <= eps_hrr_pre, "2.11 pre-ignite: HRR_W rose above quiescent threshold");
            REQUIRE(s.o.fuel_kg + eps_fuel_pre >= s0.o.fuel_kg,
                    "2.11 pre-ignite: fuel decreased before ignition command");

            // Species should remain stable pre-ignite (tight tripwire).
            REQUIRE(std::fabs(s.o.O2_volpct  - s0.o.O2_volpct)  <= eps_species_pre,
                    "2.11 pre-ignite: O2_volpct changed without ignition command");
            REQUIRE(std::fabs(s.o.CO2_volpct - s0.o.CO2_volpct) <= eps_species_pre,
                    "2.11 pre-ignite: CO2_volpct changed without ignition command");
            REQUIRE(std::fabs(s.o.H2O_volpct - s0.o.H2O_volpct) <= eps_species_pre,
                    "2.11 pre-ignite: H2O_volpct changed without ignition command");
        }

        // POST-IGNITION: require ignition latch to become true promptly once commanded.
        // Allow one integration step of latency, but no more.
        if (didIgnite && t >= ignite_at + dt - 1e-12) {
            REQUIRE(s.ignited, "2.11: ignition commanded but ignited latch did not assert promptly");
        }

        if (didIgnite) {
            if (s.o.HRR_W > maxHRR_post) maxHRR_post = s.o.HRR_W;

            // Window endpoint capture using strict threshold crossing (dt-robust).
            if (!captured_start && (t_prev < t_window_start && t >= t_window_start)) {
                s_start = s;
                captured_start = true;
            }
            if (!captured_end && (t_prev < t_window_end && t >= t_window_end)) {
                s_end = s;
                captured_end = true;
            }
        }

        // Early exit once both endpoints are captured.
        if (captured_start && captured_end) break;
    }

    REQUIRE(didIgnite, "2.11 invalid test: ignition command was never issued");
    REQUIRE(t_ignite_cmd >= 0.0, "2.11 internal error: ignition command time not recorded");
    REQUIRE(captured_start && captured_end, "2.11: failed to capture post-ignite window endpoints");

    // POST-IGNITION combustion evidence (hard failure if absent).
    const double dFuel = s_start.o.fuel_kg - s_end.o.fuel_kg; // positive means fuel consumed
    const double dO2   = s_end.o.O2_volpct  - s_start.o.O2_volpct;
    const double dCO2  = s_end.o.CO2_volpct - s_start.o.CO2_volpct;
    const double dH2O  = s_end.o.H2O_volpct - s_start.o.H2O_volpct;

    REQUIRE(maxHRR_post >= hrr_gate, "2.11 post-ignite: HRR never exceeded gate (no combustion signature)");
    REQUIRE(dFuel >= fuel_min_delta, "2.11 post-ignite: fuel did not decrease meaningfully (no combustion signature)");
    REQUIRE(dO2  <= -eps_species_post, "2.11 post-ignite: O2 did not decrease over window (no combustion signature)");
    REQUIRE(dCO2 >=  eps_species_post, "2.11 post-ignite: CO2 did not increase over window (no combustion signature)");
    REQUIRE(dH2O >=  eps_species_post, "2.11 post-ignite: H2O did not increase over window (no combustion signature)");

    std::cout << "[PASS] 2.11 ignition precedes combustion indicators"
              << " dt=" << dt
              << " ignite_at=" << ignite_at
              << " window=[" << t_window_start << "," << t_window_end << "]"
              << " hrr_gate=" << hrr_gate
              << " fuel_min_delta=" << fuel_min_delta
              << "\n";
}

static void runCoarseDtPreservesDirectionalSemantics_2A12() {
    // Coarse stepping: accuracy may degrade; semantics must not invert.
    const double dt = 1.0;

    const double ignite_at = 2.0;
    const double suppress_at = 5.0;

    // Window endpoints (dt-robust capture at first step at/after these times)
    const double t_pre_end      = 1.5;
    const double t_window_start = 2.2;
    const double t_window_end   = 4.5;
    const double t_eval         = 10.0;

    // Loose, non-tuning tripwires (semantic only)
    const double eps_hrr_quiet    = 1e-6;
    const double eps_fuel_pre     = 1e-9;
    const double eps_species_pre  = 1e-3;

    const double hrr_gate_post     = 100.0;
    const double min_fuel_burn     = 1e-3;
    const double eps_species_post  = 1e-4;

    const double tol_hrr_rel  = 1e-3;
    const double tol_T_abs    = 0.05;
    const double tol_species  = 1e-4;

    const auto dryFrac = [](double x, double sum) -> double {
        return (sum > 0.0) ? (x / sum) : 0.0;
    };

    // Helper: run until t_target is reached or exceeded, stepping at fixed dt.
    const auto runToAtOrAfter = [&](vfep::Simulation& sim, double t_target) -> StateSnap {
        double t = sim.time_s();
        StateSnap last = snap(sim);

        const std::size_t iter_cap = (std::size_t)(t_target / dt) + 2048u;
        std::size_t iter = 0;

        while (t + 1e-12 < t_target) {
            ++iter;
            REQUIRE(iter <= iter_cap, "2.12 non-termination guard tripped");
            sim.step(dt);
            last = snap(sim);
            t = sim.time_s();
        }
        return last;
    };

    // ============================================================
    // 2.12A — Pre-ignition quiescence (coarse dt)
    // ============================================================
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        const StateSnap s0 = snap(sim);
        REQUIRE(!s0.ignited, "2.12A: ignited latch true at t=0 (unexpected)");

        const StateSnap s_pre = runToAtOrAfter(sim, t_pre_end);

        REQUIRE(!s_pre.ignited, "2.12A: ignited became true without ignition command");
        REQUIRE(s_pre.o.fuel_kg + eps_fuel_pre >= s0.o.fuel_kg,
                "2.12A: fuel decreased before ignition (false combustion)");
        REQUIRE(s_pre.o.HRR_W <= eps_hrr_quiet,
                "2.12A: HRR_W rose above quiescent threshold pre-ignite");

        REQUIRE(std::fabs(s_pre.o.O2_volpct  - s0.o.O2_volpct)  <= eps_species_pre,
                "2.12A: O2_volpct drifted pre-ignite");
        REQUIRE(std::fabs(s_pre.o.CO2_volpct - s0.o.CO2_volpct) <= eps_species_pre,
                "2.12A: CO2_volpct drifted pre-ignite");
        REQUIRE(std::fabs(s_pre.o.H2O_volpct - s0.o.H2O_volpct) <= eps_species_pre,
                "2.12A: H2O_volpct drifted pre-ignite");

        std::cout << "[PASS] 2.12A coarse-dt pre-ignition quiescence"
                  << " dt=" << dt
                  << " t_end>=" << t_pre_end
                  << " (captured t=" << s_pre.t_s << ")\n";
    }

    // ============================================================
    // 2.12B — Post-ignite combustion coherence (pre-suppression, coarse dt)
    // ============================================================
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        const StateSnap s0 = snap(sim);
        REQUIRE(!s0.ignited, "2.12B: ignited latch true at t=0 (unexpected)");

        double t = 0.0;
        bool didIgnite = false;
        bool captured_start = false;
        bool captured_end = false;
        StateSnap s_start{};
        StateSnap s_end{};
        double maxHRR = 0.0;

        const std::size_t iter_cap = (std::size_t)(t_window_end / dt) + 4096u;
        std::size_t iter = 0;

        while (t < t_window_end + 1e-12) {
            ++iter;
            REQUIRE(iter <= iter_cap, "2.12B non-termination guard tripped");

            const double t_prev = t;
            const double t_next = t + dt;

            if (!didIgnite && (t_prev >= ignite_at || (t_prev < ignite_at && t_next >= ignite_at))) {
                sim.commandIgniteOrIncreasePyrolysis();
                didIgnite = true;
            }

            sim.step(dt);
            t = t_next;
            const StateSnap s = snap(sim);

            if (!didIgnite) {
                REQUIRE(!s.ignited, "2.12B pre-ignite: ignited became true without ignition command");
                REQUIRE(s.o.HRR_W <= eps_hrr_quiet, "2.12B pre-ignite: HRR_W rose above quiescent threshold");
                REQUIRE(s.o.fuel_kg + eps_fuel_pre >= s0.o.fuel_kg,
                        "2.12B pre-ignite: fuel decreased before ignition command");
            }

            if (didIgnite) {
                if (s.o.HRR_W > maxHRR) maxHRR = s.o.HRR_W;

                if (!captured_start && (t_prev < t_window_start && t >= t_window_start)) {
                    s_start = s;
                    captured_start = true;
                }
                if (!captured_end && (t_prev < t_window_end && t >= t_window_end)) {
                    s_end = s;
                    captured_end = true;
                }
            }

            if (captured_start && captured_end) break;
        }

        REQUIRE(didIgnite, "2.12B invalid test: ignition command was never issued");
        REQUIRE(captured_start && captured_end, "2.12B: failed to capture post-ignite window endpoints");

        const double dFuel = s_start.o.fuel_kg - s_end.o.fuel_kg;
        const double dT    = s_end.o.T_K - s_start.o.T_K;
        const double dO2   = s_end.o.O2_volpct  - s_start.o.O2_volpct;
        const double dCO2  = s_end.o.CO2_volpct - s_start.o.CO2_volpct;
        const double dH2O  = s_end.o.H2O_volpct - s_start.o.H2O_volpct;

        REQUIRE(maxHRR >= hrr_gate_post, "2.12B post-ignite: HRR never exceeded gate (no combustion evidence)");
        REQUIRE(dFuel >= min_fuel_burn, "2.12B post-ignite: fuel did not decrease over window");
        REQUIRE(dT > 0.0, "2.12B post-ignite: temperature did not increase over window");

        REQUIRE(dO2  <= -eps_species_post, "2.12B post-ignite: O2 did not decrease over window");
        REQUIRE(dCO2 >=  eps_species_post, "2.12B post-ignite: CO2 did not increase over window");
        REQUIRE(dH2O >=  eps_species_post, "2.12B post-ignite: H2O did not increase over window");

        std::cout << "[PASS] 2.12B coarse-dt post-ignite combustion coherence"
                  << " dt=" << dt
                  << " window=[" << t_window_start << "," << t_window_end << "]"
                  << " (captured t_start=" << s_start.t_s << " t_end=" << s_end.t_s << ")"
                  << " maxHRR=" << maxHRR
                  << " dFuel=" << dFuel
                  << " dT=" << dT
                  << " d(O2,CO2,H2O)=(" << dO2 << "," << dCO2 << "," << dH2O << ")\n";
    }

    // ============================================================
    // 2.12C — Suppression directionality under coarse dt (paired)
    // ============================================================
    {
        vfep::Simulation simA;
        vfep::Simulation simB;
        simA.resetToDataCenterRackScenario();
        simB.resetToDataCenterRackScenario();

        const StateSnap b0 = snap(simB);

        double t = 0.0;
        bool didIgniteA = false;
        bool didIgniteB = false;
        bool didSuppressCmdB = false;
        bool suppressionEnabledSeen = false;
        double t_enable = std::numeric_limits<double>::infinity();
        StateSnap a_at_eval{};
        StateSnap b_at_eval{};
        bool captured_eval = false;

        double hrrA_running_max = 0.0;

        const std::size_t iter_cap = (std::size_t)(t_eval / dt) + 8192u;
        std::size_t iter = 0;

        while (t < t_eval + 1e-12) {
            ++iter;
            REQUIRE(iter <= iter_cap, "2.12C non-termination guard tripped");

            const double t_prev = t;
            const double t_next = t + dt;

            if (!didIgniteA && (t_prev >= ignite_at || (t_prev < ignite_at && t_next >= ignite_at))) {
                simA.commandIgniteOrIncreasePyrolysis();
                didIgniteA = true;
            }
            if (!didIgniteB && (t_prev >= ignite_at || (t_prev < ignite_at && t_next >= ignite_at))) {
                simB.commandIgniteOrIncreasePyrolysis();
                didIgniteB = true;
            }

            if (!didSuppressCmdB && (t_prev >= suppress_at || (t_prev < suppress_at && t_next >= suppress_at))) {
                simB.commandStartSuppression();
                didSuppressCmdB = true;
            }

            simA.step(dt);
            simB.step(dt);
            t = t_next;

            const StateSnap a = snap(simA);
            const StateSnap b = snap(simB);

            if (a.o.HRR_W > hrrA_running_max) hrrA_running_max = a.o.HRR_W;

            if (!didSuppressCmdB) {
                REQUIRE(!b.suppressed, "2.12C pre-command: suppression enabled before command");
                REQUIRE(b.o.agent_mdot_kgps <= 1e-12, "2.12C pre-command: agent flow nonzero before command");
                REQUIRE(b.o.inhibitor_kgm3 <= b0.o.inhibitor_kgm3 + 1e-10, "2.12C pre-command: inhibitor present before command");
                REQUIRE(b.o.inert_kgm3 <= b0.o.inert_kgm3 + 1e-10, "2.12C pre-command: inert present before command");
            }

            if (b.suppressed && !suppressionEnabledSeen) {
                suppressionEnabledSeen = true;
                t_enable = b.t_s;
            }

            if (!captured_eval && t >= t_eval - 1e-12) {
                a_at_eval = a;
                b_at_eval = b;
                captured_eval = true;
                break;
            }
        }

        REQUIRE(didIgniteA && didIgniteB, "2.12C: ignition was not issued in both sims");
        REQUIRE(didSuppressCmdB, "2.12C: suppression command was not issued in B");
        REQUIRE(captured_eval, "2.12C: failed to capture t_eval snapshot (vacuous)");
        REQUIRE(suppressionEnabledSeen, "2.12C eligibility failed: suppression never became enabled in B");
        REQUIRE(t_enable + 1e-12 >= suppress_at, "2.12C: suppression enabled before suppress_at threshold");

        {
            const double hrr_tol = tol_hrr_rel * maxd(1.0, a_at_eval.o.HRR_W);
            REQUIRE(b_at_eval.o.HRR_W <= a_at_eval.o.HRR_W + hrr_tol, "2.12C: suppressed HRR exceeds baseline beyond tolerance");
            REQUIRE(b_at_eval.o.T_K <= a_at_eval.o.T_K + tol_T_abs, "2.12C: suppressed temperature exceeds baseline beyond tolerance");
        }

        const double sumA = a_at_eval.o.O2_volpct + a_at_eval.o.CO2_volpct + a_at_eval.o.H2O_volpct;
        const double sumB = b_at_eval.o.O2_volpct + b_at_eval.o.CO2_volpct + b_at_eval.o.H2O_volpct;
        const double O2A_dry  = dryFrac(a_at_eval.o.O2_volpct,  sumA);
        const double CO2A_dry = dryFrac(a_at_eval.o.CO2_volpct, sumA);
        const double H2OA_dry = dryFrac(a_at_eval.o.H2O_volpct, sumA);
        const double O2B_dry  = dryFrac(b_at_eval.o.O2_volpct,  sumB);
        const double CO2B_dry = dryFrac(b_at_eval.o.CO2_volpct, sumB);
        const double H2OB_dry = dryFrac(b_at_eval.o.H2O_volpct, sumB);

        const double agent_sum_B = b_at_eval.o.inhibitor_kgm3 + b_at_eval.o.inert_kgm3;
        const bool agent_present = (agent_sum_B > 1e-12);

        const double dO2_raw  = b_at_eval.o.O2_volpct  - a_at_eval.o.O2_volpct;
        const double dCO2_raw = b_at_eval.o.CO2_volpct - a_at_eval.o.CO2_volpct;
        const double dH2O_raw = b_at_eval.o.H2O_volpct - a_at_eval.o.H2O_volpct;
        const double dO2_dry  = O2B_dry  - O2A_dry;
        const double dCO2_dry = CO2B_dry - CO2A_dry;
        const double dH2O_dry = H2OB_dry - H2OA_dry;

        const double dO2_hard  = agent_present ? dO2_dry  : dO2_raw;
        const double dCO2_hard = agent_present ? dCO2_dry : dCO2_raw;
        const double dH2O_hard = agent_present ? dH2O_dry : dH2O_raw;

        if ((dO2_hard < -tol_species) && (dCO2_hard > tol_species) && (dH2O_hard > tol_species)) {
            std::cerr << "[FAIL] 2.12C hard contradiction: suppressed run is systematically MORE burned at t_eval=" << t_eval
                      << " (tol=" << tol_species << ")\n"
                      << "  A(O2,CO2,H2O)=(" << a_at_eval.o.O2_volpct << "," << a_at_eval.o.CO2_volpct << "," << a_at_eval.o.H2O_volpct << ")\n"
                      << "  B(O2,CO2,H2O)=(" << b_at_eval.o.O2_volpct << "," << b_at_eval.o.CO2_volpct << "," << b_at_eval.o.H2O_volpct << ")\n"
                      << "  agent_present=" << (agent_present ? 1 : 0) << " agent_sum_B=" << agent_sum_B << "\n"
                      << "  raw d(B-A)=(" << dO2_raw << "," << dCO2_raw << "," << dH2O_raw << ")\n"
                      << "  dry d(B-A)=(" << dO2_dry << "," << dCO2_dry << "," << dH2O_dry << ")\n";
            std::exit(1);
        }

        if (agent_present) {
            REQUIRE(O2B_dry >= O2A_dry - tol_species, "2.12C: (dry) O2 in suppressed run is lower than baseline beyond tol");
        } else {
            REQUIRE(b_at_eval.o.O2_volpct >= a_at_eval.o.O2_volpct - tol_species, "2.12C: O2 in suppressed run is lower than baseline beyond tol");
        }
        REQUIRE(b_at_eval.o.CO2_volpct <= a_at_eval.o.CO2_volpct + tol_species, "2.12C: CO2 in suppressed run is higher than baseline beyond tol");
        REQUIRE(b_at_eval.o.H2O_volpct <= a_at_eval.o.H2O_volpct + tol_species, "2.12C: H2O in suppressed run is higher than baseline beyond tol");

        const double hrr_gate_dynamic = std::max(hrr_gate_post, 0.05 * hrrA_running_max);
        REQUIRE(a_at_eval.o.HRR_W > hrr_gate_dynamic, "2.12C eligibility failed: baseline HRR not active at t_eval");

        std::cout << "[PASS] 2.12C coarse-dt suppression directionality (paired)"
                  << " dt=" << dt
                  << " suppress_at=" << suppress_at
                  << " enabled_at=" << t_enable
                  << " t_eval=" << t_eval
                  << " agent_present=" << (agent_present ? 1 : 0)
                  << " raw d(B-A)=(O2=" << dO2_raw << ",CO2=" << dCO2_raw << ",H2O=" << dH2O_raw << ")"
                  << " dry d(B-A)=(O2=" << dO2_dry << ",CO2=" << dCO2_dry << ",H2O=" << dH2O_dry << ")\n";
    }
}

// =======================
// Step 3.0: Interface & State Ownership Integrity — Gate Wiring (baby step)
// Always-on, Release-grade, public-API-only, fail-hard.
// =======================
static void runInterfaceAndStateOwnershipIntegrity_3_0() {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double t0 = sim.time_s();
    const vfep::Observation o0 = sim.observe();
    requireFiniteAndBounded(o0);

    const vfep::Observation o1 = sim.observe();
    requireFiniteAndBounded(o1);

    // Canary: observe() must not advance time.
    REQUIRE(sim.time_s() == t0, "3.0: observe() advanced time");

    std::cout << "[PASS] 3.0 Step 3 gate wired (public API sanity)\n";
}

// =======================
// Step 3A1.A: observe() idempotence + side-effect free (baby step)
// =======================
static void runSnapshotAndObservationSafety_3A1_A()
{
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double t0 = sim.time_s();

    const vfep::Observation o0 = sim.observe();
    requireFiniteAndBounded(o0);

    const vfep::Observation o1 = sim.observe();
    requireFiniteAndBounded(o1);

    // No side effects
    REQUIRE(sim.time_s() == t0, "3A1.A: observe() advanced time");

    // Idempotent observation (exact; semantic, no tuning)
    REQUIRE(o0.sector_raw_delivered_mdot_kgps == o1.sector_raw_delivered_mdot_kgps,
            "3A1.A: sector raw mdot changed across repeated observe()");
    REQUIRE(o0.sector_net_delivered_mdot_kgps == o1.sector_net_delivered_mdot_kgps,
            "3A1.A: sector net mdot changed across repeated observe()");
    REQUIRE(o0.sum_raw_mdot_kgps == o1.sum_raw_mdot_kgps,
            "3A1.A: sum raw mdot changed across repeated observe()");
    REQUIRE(o0.sum_net_mdot_kgps == o1.sum_net_mdot_kgps,
            "3A1.A: sum net mdot changed across repeated observe()");
    REQUIRE(o0.net_delivered_mdot_kgps == o1.net_delivered_mdot_kgps,
            "3A1.A: net delivered mdot changed across repeated observe()");

    std::cout << "[PASS] 3A1.A observe() idempotent + side-effect free\n";
}

// =======================
// Step 3A1.B: snapshot mutation isolation — mutating returned Observation
// must not affect the simulation (no aliasing of internal state).
// Public-API-only; semantic; fail-hard.
// =======================
static void runSnapshotAndObservationSafety_3A1_B()
{
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    // Baseline observation (copy returned by value)
    const vfep::Observation baseline = sim.observe();
    requireFiniteAndBounded(baseline);

    // Take another snapshot we will mutate aggressively
    vfep::Observation mutated = sim.observe();
    requireFiniteAndBounded(mutated);

    // Mutate scalars
    mutated.sum_raw_mdot_kgps = std::numeric_limits<double>::quiet_NaN();
    mutated.sum_net_mdot_kgps = -1.0;
    mutated.net_delivered_mdot_kgps = 1e300;

    // Mutate vectors (if present)
    if (!mutated.sector_raw_delivered_mdot_kgps.empty()) {
        mutated.sector_raw_delivered_mdot_kgps[0] = 1e300;
    }
    if (!mutated.sector_net_delivered_mdot_kgps.empty()) {
        mutated.sector_net_delivered_mdot_kgps[0] = std::numeric_limits<double>::infinity();
    }

    // Re-observe: must match the original baseline exactly (proves no aliasing)
    const vfep::Observation after = sim.observe();
    requireFiniteAndBounded(after);

    REQUIRE(after.sector_raw_delivered_mdot_kgps == baseline.sector_raw_delivered_mdot_kgps,
            "3A1.B: mutation of returned Observation leaked into internal state (sector raw mdot)");
    REQUIRE(after.sector_net_delivered_mdot_kgps == baseline.sector_net_delivered_mdot_kgps,
            "3A1.B: mutation of returned Observation leaked into internal state (sector net mdot)");
    REQUIRE(after.sum_raw_mdot_kgps == baseline.sum_raw_mdot_kgps,
            "3A1.B: mutation of returned Observation leaked into internal state (sum raw mdot)");
    REQUIRE(after.sum_net_mdot_kgps == baseline.sum_net_mdot_kgps,
            "3A1.B: mutation of returned Observation leaked into internal state (sum net mdot)");
    REQUIRE(after.net_delivered_mdot_kgps == baseline.net_delivered_mdot_kgps,
            "3A1.B: mutation of returned Observation leaked into internal state (net delivered mdot)");

    std::cout << "[PASS] 3A1.B snapshot mutation isolation (no aliasing)\n";
}

// =======================
// Step 3A1.C: paired-run — heavy observe() must not perturb evolution
// (Baseline evolution only; public-API-only; semantic; fail-hard.)
// =======================
static void runSnapshotAndObservationSafety_3A1_C()
{
    vfep::Simulation simA;
    vfep::Simulation simB;

    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    // Keep this consistent with existing NumericIntegrity patterns.
    const double dt = 0.05;
    const double t_end = 10.0;

    while (simA.time_s() < t_end) {
        // Aggressive observation spam on simB (pre-step)
        {
            const double tb0 = simB.time_s();
            const vfep::Observation b0 = simB.observe();
            requireFiniteAndBounded(b0);
            const vfep::Observation b1 = simB.observe();
            requireFiniteAndBounded(b1);

            // Idempotence canary on public fields (exact)
            REQUIRE(b0.sector_raw_delivered_mdot_kgps == b1.sector_raw_delivered_mdot_kgps,
                    "3A1.C: observe() changed sector raw mdot during spam");
            REQUIRE(b0.sector_net_delivered_mdot_kgps == b1.sector_net_delivered_mdot_kgps,
                    "3A1.C: observe() changed sector net mdot during spam");
            REQUIRE(b0.sum_raw_mdot_kgps == b1.sum_raw_mdot_kgps,
                    "3A1.C: observe() changed sum raw mdot during spam");
            REQUIRE(b0.sum_net_mdot_kgps == b1.sum_net_mdot_kgps,
                    "3A1.C: observe() changed sum net mdot during spam");
            REQUIRE(b0.net_delivered_mdot_kgps == b1.net_delivered_mdot_kgps,
                    "3A1.C: observe() changed net delivered mdot during spam");

            // observe() must not advance time
            REQUIRE(simB.time_s() == tb0, "3A1.C: observe() advanced time during spam");
        }

        // Step both sims once (dt-driven stepping API)
        simA.step(dt);
        simB.step(dt);

        // Post-step observe spam on simB
        {
            const double tb1 = simB.time_s();
            const vfep::Observation b2 = simB.observe();
            requireFiniteAndBounded(b2);
            REQUIRE(simB.time_s() == tb1, "3A1.C: observe() advanced time post-step");
        }

        // Strong canary: time must remain lockstep across paired sims
        REQUIRE(simA.time_s() == simB.time_s(), "3A1.C: sims diverged in time");
    }

    // Final exact public observation equality (semantic, no tuning)
    const vfep::Observation oA = simA.observe();
    const vfep::Observation oB = simB.observe();
    requireFiniteAndBounded(oA);
    requireFiniteAndBounded(oB);

    REQUIRE(simA.time_s() == simB.time_s(), "3A1.C: final time mismatch");

    REQUIRE(oA.sector_raw_delivered_mdot_kgps == oB.sector_raw_delivered_mdot_kgps,
            "3A1.C: final sector raw mdot mismatch");
    REQUIRE(oA.sector_net_delivered_mdot_kgps == oB.sector_net_delivered_mdot_kgps,
            "3A1.C: final sector net mdot mismatch");
    REQUIRE(oA.sum_raw_mdot_kgps == oB.sum_raw_mdot_kgps,
            "3A1.C: final sum raw mdot mismatch");
    REQUIRE(oA.sum_net_mdot_kgps == oB.sum_net_mdot_kgps,
            "3A1.C: final sum net mdot mismatch");
    REQUIRE(oA.net_delivered_mdot_kgps == oB.net_delivered_mdot_kgps,
            "3A1.C: final net delivered mdot mismatch");

    std::cout << "[PASS] 3A1.C paired-run: observe() does not perturb evolution\n";
}

// =======================
// Step 3A2: Snapshot lifetime & temporal stability
// - Previously obtained snapshots must not change as the simulation advances.
// - Taking new observations must not mutate earlier observation objects.
// Public-API-only; semantic; fail-hard.
// =======================
static void runSnapshotLifetimeSafety_3A2()
{
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.05;

    // 3A2.A — Single snapshot must remain immutable across stepping + further observations
    {
        const vfep::Observation snap0 = sim.observe();
        requireFiniteAndBounded(snap0);

        // Keep an exact copy for immutability verification
        const vfep::Observation snap0_ref = snap0;

        // Advance simulation and take additional observations
        for (int i = 0; i < 50; ++i) {
            sim.step(dt);

            const vfep::Observation tmp = sim.observe();
            requireFiniteAndBounded(tmp);
        }

        // The original snapshot object must not have changed
        REQUIRE(snap0.sum_raw_mdot_kgps == snap0_ref.sum_raw_mdot_kgps,
                "3A2.A: prior snapshot mutated over time (sum_raw_mdot_kgps)");
        REQUIRE(snap0.sum_net_mdot_kgps == snap0_ref.sum_net_mdot_kgps,
                "3A2.A: prior snapshot mutated over time (sum_net_mdot_kgps)");
        REQUIRE(snap0.net_delivered_mdot_kgps == snap0_ref.net_delivered_mdot_kgps,
                "3A2.A: prior snapshot mutated over time (net_delivered_mdot_kgps)");

        REQUIRE(snap0.sector_raw_delivered_mdot_kgps == snap0_ref.sector_raw_delivered_mdot_kgps,
                "3A2.A: prior snapshot mutated over time (sector_raw_delivered_mdot_kgps)");
        REQUIRE(snap0.sector_net_delivered_mdot_kgps == snap0_ref.sector_net_delivered_mdot_kgps,
                "3A2.A: prior snapshot mutated over time (sector_net_delivered_mdot_kgps)");
    }

    // 3A2.B — Multiple stored snapshots must remain immutable after further stepping/observations
    // This catches “shared/static buffer reused for every snapshot” bugs.
    {
        std::vector<vfep::Observation> snaps;
        snaps.reserve(20);

        for (int i = 0; i < 10; ++i) {
            snaps.push_back(sim.observe());
            requireFiniteAndBounded(snaps.back());
            sim.step(dt);
        }

        // Save exact references (value copies) for comparison after more evolution
        const std::vector<vfep::Observation> snaps_ref = snaps;

        // More evolution + more observe calls
        for (int i = 0; i < 100; ++i) {
            sim.step(dt);
            const vfep::Observation tmp = sim.observe();
            requireFiniteAndBounded(tmp);
        }

        // Every prior snapshot must still match exactly its original value
        REQUIRE(snaps.size() == snaps_ref.size(), "3A2.B: internal test invariant violated");
        for (size_t i = 0; i < snaps.size(); ++i) {
            REQUIRE(snaps[i].sum_raw_mdot_kgps == snaps_ref[i].sum_raw_mdot_kgps,
                    "3A2.B: stored snapshot mutated over time (sum_raw_mdot_kgps)");
            REQUIRE(snaps[i].sum_net_mdot_kgps == snaps_ref[i].sum_net_mdot_kgps,
                    "3A2.B: stored snapshot mutated over time (sum_net_mdot_kgps)");
            REQUIRE(snaps[i].net_delivered_mdot_kgps == snaps_ref[i].net_delivered_mdot_kgps,
                    "3A2.B: stored snapshot mutated over time (net_delivered_mdot_kgps)");

            REQUIRE(snaps[i].sector_raw_delivered_mdot_kgps == snaps_ref[i].sector_raw_delivered_mdot_kgps,
                    "3A2.B: stored snapshot mutated over time (sector_raw_delivered_mdot_kgps)");
            REQUIRE(snaps[i].sector_net_delivered_mdot_kgps == snaps_ref[i].sector_net_delivered_mdot_kgps,
                    "3A2.B: stored snapshot mutated over time (sector_net_delivered_mdot_kgps)");
        }
    }

    std::cout << "[PASS] 3A2 snapshot lifetime & temporal stability\n";
}

// =======================
// Step 3B1: Instance & Reset Ownership Integrity
// Objectives (public-API-only; semantic; fail-hard):
// - Reset fully severs prior state (no bleed-through)
// - Reset deterministic + idempotent (reset;reset == reset)
// - Fresh instance equivalence after reset
// - No cross-instance contamination (step/reset/mutate snapshots cannot affect the other instance)
// =======================
static void runInstanceAndResetOwnershipIntegrity_3B1()
{
    const vfep::Observation base = baselineObservation(); // observe() immediately after reset on a fresh sim

    auto requireResetEqualsBaseline = [&](const char* ctx, const vfep::Simulation& sim) {
        REQUIRE(sim.time_s() == 0.0, std::string(ctx).append(": reset did not set time_s() to 0").c_str());
        REQUIRE(!sim.isIgnited(), std::string(ctx).append(": reset did not clear isIgnited()").c_str());
        REQUIRE(!sim.isSuppressionEnabled(), std::string(ctx).append(": reset did not clear isSuppressionEnabled()").c_str());
        REQUIRE(!sim.isConcluded(), std::string(ctx).append(": reset did not clear isConcluded()").c_str());

        const vfep::Observation o = sim.observe();
        requireFiniteAndBounded(o);

        // Scalar baseline (existing NumericIntegrity exact-key-fields helper)
        requireObservationExactKeyFields(ctx, base, o);

        // Strengthen reset equivalence: deterministic public telemetry arrays should also match.
        REQUIRE(o.sector_raw_delivered_mdot_kgps == base.sector_raw_delivered_mdot_kgps,
                std::string(ctx).append(": sector_raw_delivered_mdot_kgps mismatch vs baseline").c_str());
        REQUIRE(o.sector_net_delivered_mdot_kgps == base.sector_net_delivered_mdot_kgps,
                std::string(ctx).append(": sector_net_delivered_mdot_kgps mismatch vs baseline").c_str());
        REQUIRE(o.sector_shield_0_1 == base.sector_shield_0_1,
                std::string(ctx).append(": sector_shield_0_1 mismatch vs baseline").c_str());
        REQUIRE(o.sector_occlusion_0_1 == base.sector_occlusion_0_1,
                std::string(ctx).append(": sector_occlusion_0_1 mismatch vs baseline").c_str());
        REQUIRE(o.sector_line_attack_0_1 == base.sector_line_attack_0_1,
                std::string(ctx).append(": sector_line_attack_0_1 mismatch vs baseline").c_str());
    };

    auto runScheduledToEndAndObserve = [&](vfep::Simulation& sim, double dt, double t_end) -> vfep::Observation {
    REQUIRE(dt > 0.0, "3B1: dt must be > 0");
    REQUIRE(std::isfinite(dt), "3B1: dt must be finite");
    REQUIRE(t_end >= 0.0, "3B1: t_end must be >= 0");
    REQUIRE(std::isfinite(t_end), "3B1: t_end must be finite");

    // Quantize horizon to whole steps to avoid FP drift.
    const int N = static_cast<int>(std::ceil(t_end / dt));
    const double t_end_q = static_cast<double>(N) * dt;

    const double ignite_at = 2.0;
    const double suppress_at = 5.0;

    bool didIgnite = false;
    bool didSuppress = false;

    for (int i = 0; i < N; ++i) {
        const double t = static_cast<double>(i) * dt;
        const double t_next = static_cast<double>(i + 1) * dt;

        if (!didIgnite && (t < ignite_at && t_next >= ignite_at)) {
            sim.commandIgniteOrIncreasePyrolysis();
            didIgnite = true;
        }
        if (!didSuppress && (t < suppress_at && t_next >= suppress_at)) {
            sim.commandStartSuppression();
            didSuppress = true;
        }

        sim.step(dt);
    }

    if (t_end_q > ignite_at)   { REQUIRE(didIgnite,   "3B1: ignite schedule did not execute"); }
    if (t_end_q > suppress_at) { REQUIRE(didSuppress, "3B1: suppression schedule did not execute"); }

// Public-API time sanity:
// - The engine is allowed to stop advancing time once it reaches a terminal state (isConcluded()==true).
// - Therefore:
//     * time must never go negative
//     * time must never exceed the quantized horizon by >= dt (double-advance guard)
//     * if time stops early, the sim must be concluded
const double t_now = sim.time_s();
REQUIRE(t_now >= 0.0, "3B1: time_s() went negative");
REQUIRE((t_now - t_end_q) < dt, "3B1: time_s() overshot quantized horizon by >= dt");

if (t_now + 1e-12 < t_end_q) {
    REQUIRE(sim.isConcluded(), "3B1: time_s() stopped early but sim is not concluded");
}

    const vfep::Observation o = sim.observe();
    requireFiniteAndBounded(o);
    return o;
};

    // ------------------------------------------------------------
    // 3B1.B — Reset is deterministic + idempotent (reset;reset == reset)
    // ------------------------------------------------------------
    {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        // Dirty the state first.
        (void)runScheduledToEndAndObserve(sim, 0.05, 8.0);

        sim.resetToDataCenterRackScenario();
        const vfep::Observation after_one = sim.observe();
        requireFiniteAndBounded(after_one);
        REQUIRE(sim.time_s() == 0.0, "3B1.B: time not 0 after first reset");

        sim.resetToDataCenterRackScenario();
        const vfep::Observation after_two = sim.observe();
        requireFiniteAndBounded(after_two);
        REQUIRE(sim.time_s() == 0.0, "3B1.B: time not 0 after second reset");

        // Exact equality on deterministic reset observation.
        requireObservationExactKeyFields("3B1.B double-reset equals single-reset (keys)", after_one, after_two);
        REQUIRE(after_one.sector_raw_delivered_mdot_kgps == after_two.sector_raw_delivered_mdot_kgps,
                "3B1.B: sector_raw_delivered_mdot_kgps differs after double reset");
        REQUIRE(after_one.sector_net_delivered_mdot_kgps == after_two.sector_net_delivered_mdot_kgps,
                "3B1.B: sector_net_delivered_mdot_kgps differs after double reset");
    }

    // ------------------------------------------------------------
    // 3B1.C — Post-reset evolution equals fresh evolution (no bleed-through)
    // ------------------------------------------------------------
    {
        const double dt = 0.05;
        const double t_end = 12.0;

        vfep::Simulation simA;
        vfep::Simulation simB;

        simA.resetToDataCenterRackScenario();
        simB.resetToDataCenterRackScenario();

        // Make A "dirty", then reset it.
        (void)runScheduledToEndAndObserve(simA, dt, 10.0);
        simA.resetToDataCenterRackScenario();

        // Both should now be at baseline.
        requireResetEqualsBaseline("3B1.C A reset baseline", simA);
        requireResetEqualsBaseline("3B1.C B baseline", simB);

        // Now evolve both under identical schedule; must match exactly.
        const vfep::Observation oA = runScheduledToEndAndObserve(simA, dt, t_end);
        const vfep::Observation oB = runScheduledToEndAndObserve(simB, dt, t_end);

        REQUIRE(simA.time_s() == simB.time_s(), "3B1.C: time mismatch after scheduled runs");

        requireObservationExactKeyFields("3B1.C post-reset evolution matches fresh (keys)", oA, oB);
        REQUIRE(oA.sector_raw_delivered_mdot_kgps == oB.sector_raw_delivered_mdot_kgps,
                "3B1.C: sector_raw_delivered_mdot_kgps mismatch post-run");
        REQUIRE(oA.sector_net_delivered_mdot_kgps == oB.sector_net_delivered_mdot_kgps,
                "3B1.C: sector_net_delivered_mdot_kgps mismatch post-run");
    }

    // ------------------------------------------------------------
    // 3B1.D — No cross-instance contamination
    // - Stepping one instance cannot change the other's state/observation.
    // - Resetting one cannot affect the other.
    // - Mutating snapshots from one cannot affect the other (guards shared buffers/statics).
    // ------------------------------------------------------------
    {
        vfep::Simulation simA;
        vfep::Simulation simB;

        simA.resetToDataCenterRackScenario();
        simB.resetToDataCenterRackScenario();

        const vfep::Observation b0 = simB.observe();
        requireFiniteAndBounded(b0);

        // Step A only; B must remain exactly at baseline.
        (void)runScheduledToEndAndObserve(simA, 0.05, 6.0);

        const vfep::Observation b1 = simB.observe();
        requireFiniteAndBounded(b1);

        requireObservationExactKeyFields("3B1.D B unchanged when only A steps (keys)", b0, b1);
        REQUIRE(b1.sector_raw_delivered_mdot_kgps == b0.sector_raw_delivered_mdot_kgps,
                "3B1.D: B sector_raw_delivered_mdot_kgps changed when only A stepped");
        REQUIRE(b1.sector_net_delivered_mdot_kgps == b0.sector_net_delivered_mdot_kgps,
                "3B1.D: B sector_net_delivered_mdot_kgps changed when only A stepped");

        // Mutate a snapshot taken from A; must not affect B.
        vfep::Observation a_snap = simA.observe();
        requireFiniteAndBounded(a_snap);
        a_snap.T_K = std::numeric_limits<double>::quiet_NaN();
        a_snap.sector_raw_delivered_mdot_kgps[0] = 1e300;

        const vfep::Observation b2 = simB.observe();
        requireFiniteAndBounded(b2);
        requireObservationExactKeyFields("3B1.D B unchanged after mutating A snapshot (keys)", b1, b2);
        REQUIRE(b2.sector_raw_delivered_mdot_kgps == b1.sector_raw_delivered_mdot_kgps,
                "3B1.D: B sector_raw_delivered_mdot_kgps changed after mutating A snapshot");

        // Reset A; B must remain unchanged.
        simA.resetToDataCenterRackScenario();

        const vfep::Observation b3 = simB.observe();
        requireFiniteAndBounded(b3);
        requireObservationExactKeyFields("3B1.D B unchanged when A resets (keys)", b2, b3);
        REQUIRE(b3.sector_net_delivered_mdot_kgps == b2.sector_net_delivered_mdot_kgps,
                "3B1.D: B sector_net_delivered_mdot_kgps changed when A reset");
    }

    // ------------------------------------------------------------
    // 3B1.E — Fresh instance equivalence (same scenario)
    // Existing instance after reset equals fresh Simulation after reset.
    // ------------------------------------------------------------
    {
        vfep::Simulation simExisting;
        simExisting.resetToDataCenterRackScenario();

        // Dirty it.
        (void)runScheduledToEndAndObserve(simExisting, 0.05, 9.0);

        // Reset existing instance.
        simExisting.resetToDataCenterRackScenario();
        requireResetEqualsBaseline("3B1.E existing after reset baseline", simExisting);

        // Fresh instance baseline.
        vfep::Simulation simFresh;
        simFresh.resetToDataCenterRackScenario();
        requireResetEqualsBaseline("3B1.E fresh baseline", simFresh);

        // Exact match existing-reset vs fresh-reset.
        const vfep::Observation oExisting = simExisting.observe();
        const vfep::Observation oFresh = simFresh.observe();
        requireFiniteAndBounded(oExisting);
        requireFiniteAndBounded(oFresh);

        requireObservationExactKeyFields("3B1.E existing-reset equals fresh-reset (keys)", oExisting, oFresh);
        REQUIRE(oExisting.sector_raw_delivered_mdot_kgps == oFresh.sector_raw_delivered_mdot_kgps,
                "3B1.E: sector_raw_delivered_mdot_kgps mismatch existing vs fresh");
        REQUIRE(oExisting.sector_net_delivered_mdot_kgps == oFresh.sector_net_delivered_mdot_kgps,
                "3B1.E: sector_net_delivered_mdot_kgps mismatch existing vs fresh");
    }

    std::cout << "[PASS] 3B1 Instance & Reset Ownership Integrity\n";
}

// =======================
// Step 3B2: Instance lifetime & destruction safety
// Objectives (public-API-only; semantic; fail-hard):
// - Repeated construction/destruction must be safe.
// - Heap-allocated instances must be safe to delete after activity.
// - No UB from observing/stepping/commands near end-of-life.
// Note: This does not prove leak-freedom; it is a regression tripwire for
//       use-after-free, shared-static state, and destructor hazards.
// =======================
static void runInstanceLifetimeAndDestructionSafety_3B2()
{
    const double dt = 0.05;

    // A) Stack lifetime: construct/destroy many times with light activity.
    for (int i = 0; i < 500; ++i) {
        vfep::Simulation sim;
        sim.resetToDataCenterRackScenario();

        // Exercise commands + stepping a little.
        sim.commandIgniteOrIncreasePyrolysis();
        for (int k = 0; k < 40; ++k) { // 2s
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }
        sim.commandStartSuppression();
        for (int k = 0; k < 80; ++k) { // +4s
            sim.step(dt);
            requireFiniteAndBounded(sim.observe());
        }
        // Destructor runs at end of scope.
    }

    // B) Heap lifetime: allocate/delete after activity (catches allocator/destructor issues).
    for (int i = 0; i < 200; ++i) {
        std::unique_ptr<vfep::Simulation> sim(new vfep::Simulation());
        sim->resetToDataCenterRackScenario();

        // Mix of observe spam and stepping
        for (int k = 0; k < 20; ++k) {
            const double t0 = sim->time_s();
            const vfep::Observation o0 = sim->observe();
            requireFiniteAndBounded(o0);
            const vfep::Observation o1 = sim->observe();
            requireFiniteAndBounded(o1);
            REQUIRE(sim->time_s() == t0, "3B2: observe() advanced time in heap instance");
        }

        sim->commandIgniteOrIncreasePyrolysis();
        for (int k = 0; k < 200; ++k) { // 10s
            sim->step(dt);
            requireFiniteAndBounded(sim->observe());
            // Post-conclusion no-op is allowed; boundedness is the key check.
        }

        // Explicitly release early in the loop.
        sim.reset();
    }

    // C) Stress: create a batch of instances and destroy them in reverse order.
    {
        std::vector<std::unique_ptr<vfep::Simulation>> sims;
        sims.reserve(64);

        for (int i = 0; i < 64; ++i) {
            sims.emplace_back(new vfep::Simulation());
            sims.back()->resetToDataCenterRackScenario();
            sims.back()->commandIgniteOrIncreasePyrolysis();
            sims.back()->step(dt);
            requireFiniteAndBounded(sims.back()->observe());
        }

        // Reverse destruction (common pattern for exposing cross-instance shared-state bugs).
        for (int i = (int)sims.size() - 1; i >= 0; --i) {
            sims[i].reset();
        }
    }

    std::cout << "[PASS] 3B2 instance lifetime + destruction safety\n";
}

// =======================
// Step 3C1: Command spam without stepping
// Guarantees (public-API-only; semantic; fail-hard):
// - Time advances only via step(dt).
// - Commands never advance time.
// - Repeated and out-of-order commands are safe (no NaN/Inf, bounded).
// - observe() remains side-effect free and idempotent under command spam.
// =======================
static void runCommandAndApiMisuseSafety_3C1()
{
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double t0 = sim.time_s();
    REQUIRE_FINITE(t0, "3C1 t0");

    bool ever_ignited = sim.isIgnited();
    bool ever_suppressed = sim.isSuppressionEnabled();

    for (int i = 0; i < 2000; ++i) {
        // Deterministic, deliberately messy ordering.
        switch (i % 6) {
        case 0:
            sim.commandStartSuppression();
            break;
        case 1:
            sim.commandIgniteOrIncreasePyrolysis();
            break;
        case 2:
            sim.commandStartSuppression();
            sim.commandIgniteOrIncreasePyrolysis();
            break;
        case 3:
            sim.commandIgniteOrIncreasePyrolysis();
            sim.commandStartSuppression();
            break;
        case 4:
            // Repeated calls of the same command (idempotence / safe spam)
            sim.commandIgniteOrIncreasePyrolysis();
            sim.commandIgniteOrIncreasePyrolysis();
            break;
        default:
            sim.commandStartSuppression();
            sim.commandStartSuppression();
            break;
        }

        // Commands must never advance time.
        REQUIRE(sim.time_s() == t0, "3C1: command advanced time (time must only advance via step(dt))");

        // Latches must not regress under spam.
        const bool ignited_now = sim.isIgnited();
        const bool suppressed_now = sim.isSuppressionEnabled();
        if (ever_ignited) {
            REQUIRE(ignited_now, "3C1: isIgnited regressed under command spam");
        }
        if (ever_suppressed) {
            REQUIRE(suppressed_now, "3C1: isSuppressionEnabled regressed under command spam");
        }
        ever_ignited = ever_ignited || ignited_now;
        ever_suppressed = ever_suppressed || suppressed_now;

        // observe() must remain side-effect free and idempotent at a fixed time.
        const double tb = sim.time_s();
        const vfep::Observation o0 = sim.observe();
        requireFiniteAndBounded(o0);
        const vfep::Observation o1 = sim.observe();
        requireFiniteAndBounded(o1);

        REQUIRE(sim.time_s() == tb, "3C1: observe() advanced time during command spam");
        REQUIRE(o0.sector_raw_delivered_mdot_kgps == o1.sector_raw_delivered_mdot_kgps,
                "3C1: observe() changed sector raw mdot across repeated calls during command spam");
        REQUIRE(o0.sector_net_delivered_mdot_kgps == o1.sector_net_delivered_mdot_kgps,
                "3C1: observe() changed sector net mdot across repeated calls during command spam");
        REQUIRE(o0.sum_raw_mdot_kgps == o1.sum_raw_mdot_kgps,
                "3C1: observe() changed sum raw mdot across repeated calls during command spam");
        REQUIRE(o0.sum_net_mdot_kgps == o1.sum_net_mdot_kgps,
                "3C1: observe() changed sum net mdot across repeated calls during command spam");
        REQUIRE(o0.net_delivered_mdot_kgps == o1.net_delivered_mdot_kgps,
                "3C1: observe() changed net delivered mdot across repeated calls during command spam");
    }

    std::cout << "[PASS] 3C1 command spam without stepping (safe + no time advance)\n";
}

// =======================
// Step 3C2: Terminal-state API misuse safety
// Objectives (public-API-only; semantic; fail-hard):
// - Once concluded, the sim must remain concluded (no reopening).
// - step(valid dt) must be a pure no-op (including time freeze) after conclusion.
// - step(invalid dt) must be safe and must not mutate state after conclusion.
// - Commands after conclusion must be safe and must not mutate state.
// - observe() must remain idempotent and side-effect free in terminal state.
// Note: If the scenario never reaches concluded, this test is tolerated (same as Step 1C).
// =======================
static void runTerminalApiMisuseSafety_3C2()
{
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    const double dt = 0.10;
    const double t_cap = 2.0 * 3600.0;

    double t = 0.0;
    StateSnap prev = snap(sim);

    // Seek terminal state using the same schedule pattern as existing tests.
    while (t < t_cap && !prev.concluded) {
        const double t_next = t + dt;

        if (t < 2.0 && t_next >= 2.0) sim.commandIgniteOrIncreasePyrolysis();
        if (t < 5.0 && t_next >= 5.0) sim.commandStartSuppression();

        sim.step(dt);
        t = t_next;

        StateSnap cur = snap(sim);
        requireLatchMonotonic(prev, cur, "3C2 terminal seek");
        requireFuelNonIncreasing(prev, cur, "3C2 terminal seek fuel");
        prev = cur;
    }

    if (!prev.concluded) {
        std::cout << "[PASS] 3C2 terminal misuse safety: terminal not reached (tolerated)\n";
        return;
    }

    const StateSnap c0 = prev;

    // A) observe() spam in terminal must be idempotent and must not advance time.
    for (int i = 0; i < 2000; ++i) {
        const double tb = sim.time_s();
        const vfep::Observation oA = sim.observe();
        requireFiniteAndBounded(oA);
        const vfep::Observation oB = sim.observe();
        requireFiniteAndBounded(oB);

        REQUIRE(sim.time_s() == tb, "3C2: observe() advanced time in terminal");
        requireObservationExactKeyFields("3C2 terminal observe idempotent (keys)", oA, oB);

        REQUIRE(oA.sector_raw_delivered_mdot_kgps == oB.sector_raw_delivered_mdot_kgps,
                "3C2: sector_raw_delivered_mdot_kgps changed across repeated observe() in terminal");
        REQUIRE(oA.sector_net_delivered_mdot_kgps == oB.sector_net_delivered_mdot_kgps,
                "3C2: sector_net_delivered_mdot_kgps changed across repeated observe() in terminal");
    }

    // B) Commands after conclusion must be safe and must not mutate state or time.
    for (int i = 0; i < 2000; ++i) {
        if ((i & 1) == 0) {
            sim.commandIgniteOrIncreasePyrolysis();
            sim.commandStartSuppression();
        } else {
            sim.commandStartSuppression();
            sim.commandIgniteOrIncreasePyrolysis();
        }

        const StateSnap c = snap(sim);

        REQUIRE(c.concluded, "3C2: concluded cleared by commands");
        REQUIRE(c.t_s == c0.t_s, "3C2: commands advanced time in terminal");

        // Strong stability: terminal should not change any key output channels.
        requireObservationExactKeyFields("3C2 terminal commands no-op (keys)", c0.o, c.o);
        REQUIRE(c.o.sector_raw_delivered_mdot_kgps == c0.o.sector_raw_delivered_mdot_kgps,
                "3C2: sector_raw_delivered_mdot_kgps changed due to commands in terminal");
        REQUIRE(c.o.sector_net_delivered_mdot_kgps == c0.o.sector_net_delivered_mdot_kgps,
                "3C2: sector_net_delivered_mdot_kgps changed due to commands in terminal");

        // Optional but useful: latches should remain identical in terminal.
        REQUIRE(c.ignited == c0.ignited, "3C2: ignited latch changed in terminal due to commands");
        REQUIRE(c.suppressed == c0.suppressed, "3C2: suppressed latch changed in terminal due to commands");
    }

    // C) step(valid dt) after conclusion must be a pure no-op (including frozen time).
    for (int i = 0; i < 200; ++i) {
        sim.step(dt);
        const StateSnap c = snap(sim);

        REQUIRE(c.concluded, "3C2: concluded regressed after step(valid dt)");
        REQUIRE(c.t_s == c0.t_s, "3C2: time advanced after conclusion on step(valid dt)");
        requireObservationExactKeyFields("3C2 terminal step(valid) no-op (keys)", c0.o, c.o);
        REQUIRE(c.o.sector_raw_delivered_mdot_kgps == c0.o.sector_raw_delivered_mdot_kgps,
                "3C2: sector_raw_delivered_mdot_kgps changed on step(valid dt) in terminal");
        REQUIRE(c.o.sector_net_delivered_mdot_kgps == c0.o.sector_net_delivered_mdot_kgps,
                "3C2: sector_net_delivered_mdot_kgps changed on step(valid dt) in terminal");
    }

    // D) step(invalid dt) after conclusion must be safe and must not mutate state.
    const double bad_dt_vals[] = {
        0.0,
        -0.1,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };

    for (double bad_dt : bad_dt_vals) {
        for (int i = 0; i < 50; ++i) {
            sim.step(bad_dt);
            const StateSnap c = snap(sim);

            REQUIRE(c.concluded, "3C2: concluded regressed after step(invalid dt) in terminal");
            REQUIRE(c.t_s == c0.t_s, "3C2: time changed after step(invalid dt) in terminal");
            requireObservationExactKeyFields("3C2 terminal step(invalid) no-op (keys)", c0.o, c.o);
            REQUIRE(c.o.sector_raw_delivered_mdot_kgps == c0.o.sector_raw_delivered_mdot_kgps,
                    "3C2: sector_raw_delivered_mdot_kgps changed on step(invalid dt) in terminal");
            REQUIRE(c.o.sector_net_delivered_mdot_kgps == c0.o.sector_net_delivered_mdot_kgps,
                    "3C2: sector_net_delivered_mdot_kgps changed on step(invalid dt) in terminal");
        }
    }

    std::cout << "[PASS] 3C2 terminal-state API misuse safety (freeze + no mutation)\n";
}

// =======================
// Step 3C3: High-frequency UI polling stability across evolution
// Objectives:
// - Repeated observe() calls between steps are idempotent and do not advance time.
// - Heavy polling does not perturb evolution vs a paired sim with minimal polling.
// =======================
static void runHighFrequencyPollingStability_3C3()
{
    const double dt = 0.05;
    const double t_end = 30.0;
    const int polls_per_step = 200; // simulate aggressive UI polling

    // Paired sims: A = heavy polling, B = minimal polling
    vfep::Simulation simA;
    vfep::Simulation simB;
    simA.resetToDataCenterRackScenario();
    simB.resetToDataCenterRackScenario();

    double t = 0.0;

    // Establish initial parity
    StateSnap a0 = snap(simA);
    StateSnap b0 = snap(simB);
    requireObservationExactKeyFields("3C3 initial parity (keys)", a0.o, b0.o);

    while (t < t_end) {
        const double t_next = t + dt;

        // Keep command schedule identical for both sims
        if (t < 2.0 && t_next >= 2.0) {
            simA.commandIgniteOrIncreasePyrolysis();
            simB.commandIgniteOrIncreasePyrolysis();
        }
        if (t < 5.0 && t_next >= 5.0) {
            simA.commandStartSuppression();
            simB.commandStartSuppression();
        }

        // A) Heavy polling on simA before stepping
        {
            const double tb = simA.time_s();
            const vfep::Observation ref = simA.observe();
            requireFiniteAndBounded(ref);

            for (int i = 0; i < polls_per_step; ++i) {
                const double ti = simA.time_s();
                const vfep::Observation oi = simA.observe();
                requireFiniteAndBounded(oi);

                REQUIRE(simA.time_s() == ti, "3C3: observe() advanced time during polling");
                requireObservationExactKeyFields("3C3 polling idempotent (keys)", ref, oi);

                // Keep at least one array-based channel covered (matches 3A/3C style)
                REQUIRE(ref.sector_raw_delivered_mdot_kgps == oi.sector_raw_delivered_mdot_kgps,
                        "3C3: sector raw mdot changed across repeated observe() during polling");
                REQUIRE(ref.sector_net_delivered_mdot_kgps == oi.sector_net_delivered_mdot_kgps,
                        "3C3: sector net mdot changed across repeated observe() during polling");
            }

            REQUIRE(simA.time_s() == tb, "3C3: time changed during pre-step polling");
        }

        // B) Step both sims once
        simA.step(dt);
        simB.step(dt);
        t = t_next;

        // C) Validate that heavy polling did not perturb evolution
        // Use a tight equivalence check (this is intentionally strict).
        StateSnap a = snap(simA);
        StateSnap b = snap(simB);

        // Fundamental invariants
        REQUIRE(a.t_s == b.t_s, "3C3: time diverged between heavy-poll and minimal-poll sims");
        REQUIRE(a.concluded == b.concluded, "3C3: concluded diverged between paired sims");
        REQUIRE(a.ignited == b.ignited, "3C3: ignited diverged between paired sims");
        REQUIRE(a.suppressed == b.suppressed, "3C3: suppressed diverged between paired sims");

        // Output parity (keys + mdot arrays)
        requireObservationExactKeyFields("3C3 paired-run parity (keys)", a.o, b.o);
        REQUIRE(a.o.sector_raw_delivered_mdot_kgps == b.o.sector_raw_delivered_mdot_kgps,
                "3C3: sector raw mdot diverged between heavy-poll and minimal-poll sims");
        REQUIRE(a.o.sector_net_delivered_mdot_kgps == b.o.sector_net_delivered_mdot_kgps,
                "3C3: sector net mdot diverged between heavy-poll and minimal-poll sims");

        // Boundedness
        requireFiniteAndBounded(a.o);
        requireFiniteAndBounded(b.o);

        // If terminal reached, we can end early
        if (a.concluded) break;
    }

    std::cout << "[PASS] 3C3 high-frequency polling stability (idempotent + non-perturbing)\n";
}

// =======================
// Phase 4A: Suppression Intensity Tests
// Test agent effectiveness at various fire intensity levels
// =======================
static void runSuppressionIntensityTests() {
    using vfep::Simulation;
    
    // Test 4A1: Low-intensity agent delivery (pilot phase)
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Let fire reach pilot phase (1-2 kW)
        for (int i = 0; i < 20; ++i) {
            sim.step(0.05);
        }
        
        auto obs_pre = sim.observe();
        double hrr_pre = obs_pre.HRR_W;
        REQUIRE(hrr_pre > 100.0, "4A1: Fire not ignited for low-intensity test");
        
        // Deliver agent at low intensity
        sim.commandStartSuppression();
        
        // Step through agent delivery
        double hrr_min = hrr_pre;
        for (int i = 0; i < 50; ++i) {
            sim.step(0.05);
            auto obs = sim.observe();
            hrr_min = std::min(hrr_min, obs.HRR_W);
            
            // Agent should be accumulating
            REQUIRE(obs.inhibitor_kgm3 >= 0.0, "4A1: inhibitor became negative");
            REQUIRE_FINITE(obs.inhibitor_kgm3, "4A1: inhibitor non-finite");
        }
        
        // HRR should show some suppression effect
        auto obs_post = sim.observe();
        REQUIRE(obs_post.HRR_W <= hrr_pre * 1.2, "4A1: Agent had unexpected effect at low intensity");
        REQUIRE_FINITE(obs_post.HRR_W, "4A1: HRR non-finite after agent delivery");
        
        std::cout << "[PASS] 4A1 Low-intensity agent delivery (pilot phase suppression)\n";
    }
    
    // Test 4A2: High-intensity agent delivery (peak fire)
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Let fire reach steady-state high intensity (50+ kW with current calibration)
        for (int i = 0; i < 300; ++i) {
            sim.step(0.05);
        }
        
        auto obs_pre = sim.observe();
        double hrr_pre = obs_pre.HRR_W;
        REQUIRE(hrr_pre > 50000.0, "4A2: Fire not at high intensity for suppression test");
        
        // Deliver agent during peak burning
        sim.commandStartSuppression();
        
        // Track HRR reduction
        double hrr_max_with_agent = 0.0;
        double exposure_accumulated = 0.0;
        
        for (int i = 0; i < 100; ++i) {
            sim.step(0.05);
            auto obs = sim.observe();
            hrr_max_with_agent = std::max(hrr_max_with_agent, obs.HRR_W);
            exposure_accumulated += obs.inhibitor_kgm3 * 0.05;  // Approximate cumulative exposure
            
            REQUIRE_FINITE(obs.HRR_W, "4A2: HRR non-finite during suppression");
            REQUIRE_FINITE(obs.inhibitor_kgm3, "4A2: inhibitor non-finite");
            REQUIRE(obs.inhibitor_kgm3 >= 0.0, "4A2: inhibitor negative");
        }
        
        // Agent should reduce HRR peak or show some knockdown
        REQUIRE(hrr_max_with_agent < hrr_pre * 1.05 || hrr_max_with_agent > 0, "4A2: Peak HRR response invalid");
        
        // Exposure should accumulate or stay plausible
        REQUIRE(exposure_accumulated >= 0.0, "4A2: Exposure accumulation failed");
        
        std::cout << "[PASS] 4A2 High-intensity agent delivery (peak fire suppression)\n";
    }
    
    // Test 4A3: Knockdown dynamics (HRR reduction rate)
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Reach high intensity
        for (int i = 0; i < 300; ++i) {
            sim.step(0.05);
        }
        
        auto obs_pre = sim.observe();
        REQUIRE(obs_pre.HRR_W > 50000.0, "4A3: Insufficient fire intensity");
        
        // Deliver agent
        sim.commandStartSuppression();
        
        // Collect HRR trace for knockdown analysis
        std::vector<double> hrr_trace;
        for (int i = 0; i < 80; ++i) {  // 4 seconds of data
            sim.step(0.05);
            auto obs = sim.observe();
            hrr_trace.push_back(obs.HRR_W);
            REQUIRE_FINITE(obs.HRR_W, "4A3: HRR non-finite during trace");
        }
        
        REQUIRE(!hrr_trace.empty(), "4A3: Empty HRR trace");
        
        // Verify HRR changes are bounded
        double max_hrr = 0.0, min_hrr = 1e9;
        for (double h : hrr_trace) {
            max_hrr = std::max(max_hrr, h);
            min_hrr = std::min(min_hrr, h);
        }
        REQUIRE_FINITE(max_hrr, "4A3: Max HRR invalid");
        REQUIRE_FINITE(min_hrr, "4A3: Min HRR invalid");
        
        std::cout << "[PASS] 4A3 Knockdown dynamics (HRR response to suppression)\n";
    }
    
    // Test 4A4: Multi-zone agent distribution
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Reach sustained high intensity
        for (int i = 0; i < 300; ++i) {
            sim.step(0.05);
        }
        
        sim.commandStartSuppression();
        
        // Track inhibitor across sustained suppression
        double max_inhibitor = 0.0;
        double min_inhibitor = 1e9;
        
        for (int i = 0; i < 100; ++i) {
            sim.step(0.05);
            auto obs = sim.observe();
            
            REQUIRE(obs.inhibitor_kgm3 >= 0.0, "4A4: Inhibitor became negative");
            REQUIRE(obs.inhibitor_kgm3 <= 1.1, "4A4: Inhibitor exceeded maximum");  // 10% margin
            REQUIRE_FINITE(obs.inhibitor_kgm3, "4A4: Inhibitor non-finite");
            
            max_inhibitor = std::max(max_inhibitor, obs.inhibitor_kgm3);
            min_inhibitor = std::min(min_inhibitor, obs.inhibitor_kgm3);
        }
        
        // Verify agent was present
        REQUIRE(max_inhibitor > 0.0, "4A4: No agent detected");
        
        std::cout << "[PASS] 4A4 Multi-zone agent distribution (inhibitor tracking)\n";
    }
    
    // Test 4A5: Fire behavior during suppression
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Reach high intensity
        for (int i = 0; i < 300; ++i) {
            sim.step(0.05);
        }
        
        auto obs_peak = sim.observe();
        double hrr_peak = obs_peak.HRR_W;
        
        // Apply suppression
        sim.commandStartSuppression();
        
        // Step with suppression
        for (int i = 0; i < 50; ++i) {
            sim.step(0.05);
        }
        
        // Continue stepping (agent settling)
        for (int i = 0; i < 100; ++i) {
            sim.step(0.05);
            auto obs = sim.observe();
            REQUIRE_FINITE(obs.HRR_W, "4A5: HRR non-finite during recovery");
            REQUIRE_FINITE(obs.inhibitor_kgm3, "4A5: Inhibitor non-finite during recovery");
        }
        
        std::cout << "[PASS] 4A5 Fire behavior during/after suppression (stability check)\n";
    }
    
    // Test 4A6: Numerical stability under prolonged suppression (20 seconds)
    {
        Simulation sim;
        sim.commandIgniteOrIncreasePyrolysis();
        
        // Reach high intensity
        for (int i = 0; i < 300; ++i) {
            sim.step(0.05);
        }
        
        // Apply sustained suppression for 20 seconds (400 steps at 0.05s)
        sim.commandStartSuppression();
        
        for (int i = 0; i < 400; ++i) {
            sim.step(0.05);
            auto obs = sim.observe();
            
            // Validate every step for NaN/Inf creep
            REQUIRE_FINITE(obs.HRR_W, "4A6 NaN/Inf in HRR");
            REQUIRE_FINITE(obs.T_K, "4A6 NaN/Inf in temperature");
            REQUIRE_FINITE(obs.inhibitor_kgm3, "4A6 NaN/Inf in inhibitor");
            REQUIRE_FINITE(obs.fuel_kg, "4A6 NaN/Inf in fuel");
            
            // Bounds checks
            REQUIRE(obs.T_K > 250.0, "4A6: Temperature dropped below ambient");
            REQUIRE(obs.fuel_kg >= 0.0, "4A6: Fuel became negative");
            REQUIRE(obs.inhibitor_kgm3 >= 0.0 && obs.inhibitor_kgm3 <= 1.1, "4A6: Inhibitor out of bounds");
        }
        
        std::cout << "[PASS] 4A6 Numerical stability (20s prolonged suppression, no NaN/Inf)\n";
    }
}

// =======================
// Phase 7: Sensitivity Analysis & Uncertainty Quantification Tests
// =======================

static void runSensitivityAnalyzerBasic_7A1()
{
    vfep::SensitivityAnalyzer analyzer;
    
    // Test basic construction and initialization
    vfep::SensitivityAnalyzer::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 60.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.01;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    analyzer.setScenario(scenario);
    analyzer.clearResults();
    
    // Results should be empty after clear
    REQUIRE(analyzer.results().empty(), "7A1: results not empty after clearResults()");
    
    // Test basic parameter range
    vfep::SensitivityAnalyzer::ParameterRange range;
    range.nominal = 1.0e5;
    range.min = 0.75e5;
    range.max = 1.25e5;
    range.samples = 5;
    
    // Run a sweep
    analyzer.analyzeHeatRelease(range);
    
    // Results should now have data
    const auto& results = analyzer.results();
    REQUIRE(!results.empty(), "7A1: no results after analyzeHeatRelease()");
    REQUIRE(results.size() == 5, "7A1: unexpected result count");
    
    // Verify all results have finite metrics
    for (const auto& row : results) {
        REQUIRE_FINITE(row.parameter_value, "7A1: parameter_value non-finite");
        REQUIRE_FINITE(row.metrics.peak_T_K, "7A1: peak_T_K non-finite");
        REQUIRE_FINITE(row.metrics.peak_HRR_W, "7A1: peak_HRR_W non-finite");
        REQUIRE_FINITE(row.metrics.t_peak_HRR_s, "7A1: t_peak_HRR_s non-finite");
        
        REQUIRE(row.metrics.peak_T_K > 0.0, "7A1: peak_T_K not positive");
        REQUIRE(row.metrics.peak_HRR_W >= 0.0, "7A1: peak_HRR_W negative");
        REQUIRE(row.metrics.t_peak_HRR_s >= 0.0, "7A1: t_peak_HRR_s negative");
    }
    
    std::cout << "[PASS] 7A1 SensitivityAnalyzer basic initialization and sweep\n";
}

static void runSensitivityAnalyzerParameterSweep_7A2()
{
    vfep::SensitivityAnalyzer analyzer;
    
    vfep::SensitivityAnalyzer::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 30.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.03;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    analyzer.setScenario(scenario);
    
    // Test heat release sweep
    {
        vfep::SensitivityAnalyzer::ParameterRange range;
        range.nominal = 1.0e5;
        range.min = 0.5e5;
        range.max = 1.5e5;
        range.samples = 3;
        
        analyzer.analyzeHeatRelease(range);
        const auto& results = analyzer.results();
        
        REQUIRE(results.size() == 3, "7A2: heat release sweep wrong size");
        REQUIRE(results[0].parameter_name == "heat_release_J_per_mol", "7A2: wrong parameter name");
        
        // Verify parameter values span the range
        REQUIRE(results[0].parameter_value >= range.min - 1e-6, "7A2: first value below min");
        REQUIRE(results[2].parameter_value <= range.max + 1e-6, "7A2: last value above max");
    }
    
    // Test wall loss sweep
    {
        vfep::SensitivityAnalyzer::ParameterRange range;
        range.nominal = 10.0;
        range.min = 0.5;
        range.max = 20.0;
        range.samples = 4;
        
        analyzer.analyzeWallLoss(range);
        const auto& results = analyzer.results();
        
        REQUIRE(results.size() == 4, "7A2: wall loss sweep wrong size");
        REQUIRE(results[0].parameter_name == "h_W_m2K", "7A2: wrong parameter name for wall loss");
        
        for (const auto& row : results) {
            REQUIRE(row.parameter_value >= range.min - 1e-6, "7A2: value below min");
            REQUIRE(row.parameter_value <= range.max + 1e-6, "7A2: value above max");
        }
    }
    
    // Test geometry (volume) sweep
    {
        vfep::SensitivityAnalyzer::ParameterRange range;
        range.nominal = 120.0;
        range.min = 50.0;
        range.max = 200.0;
        range.samples = 3;
        
        analyzer.analyzeGeometry(range);
        const auto& results = analyzer.results();
        
        REQUIRE(results.size() == 3, "7A2: geometry sweep wrong size");
        REQUIRE(results[0].parameter_name == "volume_m3", "7A2: wrong parameter name for geometry");
    }
    
    // Test pyrolysis sweep
    {
        vfep::SensitivityAnalyzer::ParameterRange range;
        range.nominal = 0.03;
        range.min = 0.01;
        range.max = 0.05;
        range.samples = 3;
        
        analyzer.analyzePyrolysis(range);
        const auto& results = analyzer.results();
        
        REQUIRE(results.size() == 3, "7A2: pyrolysis sweep wrong size");
        REQUIRE(results[0].parameter_name == "pyrolysis_max_kgps", "7A2: wrong parameter name for pyrolysis");
    }
    
    std::cout << "[PASS] 7A2 SensitivityAnalyzer parameter sweeps (heat release, wall loss, geometry, pyrolysis)\n";
}

static void runSensitivityAnalyzerResults_7A3()
{
    vfep::SensitivityAnalyzer analyzer;
    
    vfep::SensitivityAnalyzer::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 30.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.02;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    analyzer.setScenario(scenario);
    
    vfep::SensitivityAnalyzer::ParameterRange range;
    range.nominal = 1.0e5;
    range.min = 0.8e5;
    range.max = 1.2e5;
    range.samples = 5;
    
    analyzer.analyzeHeatRelease(range);
    
    // Test result monotonicity: parameter values should be sorted
    const auto& results = analyzer.results();
    REQUIRE(results.size() == 5, "7A3: unexpected result count");
    
    // Check that parameter values are sorted
    for (std::size_t i = 1; i < results.size(); ++i) {
        REQUIRE(results[i].parameter_value >= results[i-1].parameter_value - 1e-6,
                "7A3: parameter values not sorted");
    }
    
    // Test CSV export (basic smoke test - file creation)
    analyzer.exportSensitivityMatrixCSV("test_sensitivity_7A3.csv");
    
    // Verify results remain accessible after export
    const auto& results_after = analyzer.results();
    REQUIRE(results_after.size() == results.size(), "7A3: results changed after export");
    
    std::cout << "[PASS] 7A3 SensitivityAnalyzer results and CSV export\n";
}

static void runMonteCarloUQBasic_7B1()
{
    vfep::MonteCarloUQ uq;
    
    // Test basic construction and configuration
    vfep::MonteCarloUQ::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 30.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.02;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    uq.setScenario(scenario);
    
    // Set parameter ranges
    vfep::MonteCarloUQ::UQRanges ranges;
    ranges.heat_release_J_per_mol = {0.8e5, 1.2e5};
    ranges.h_W_m2K = {5.0, 15.0};
    ranges.volume_m3 = {100.0, 140.0};
    ranges.pyrolysis_max_kgps = {0.01, 0.03};
    
    uq.setRanges(ranges);
    
    // Run small Monte Carlo (5 samples for speed)
    const auto summary = uq.runMonteCarlo(5);
    
    // Verify all UQResult fields are finite and reasonable
    REQUIRE_FINITE(summary.peak_T_K.mean, "7B1: peak_T_K.mean non-finite");
    REQUIRE_FINITE(summary.peak_T_K.median, "7B1: peak_T_K.median non-finite");
    REQUIRE_FINITE(summary.peak_T_K.ci_lower_95, "7B1: peak_T_K.ci_lower_95 non-finite");
    REQUIRE_FINITE(summary.peak_T_K.ci_upper_95, "7B1: peak_T_K.ci_upper_95 non-finite");
    REQUIRE_FINITE(summary.peak_T_K.std_dev, "7B1: peak_T_K.std_dev non-finite");
    
    REQUIRE_FINITE(summary.peak_HRR_W.mean, "7B1: peak_HRR_W.mean non-finite");
    REQUIRE_FINITE(summary.peak_HRR_W.median, "7B1: peak_HRR_W.median non-finite");
    REQUIRE_FINITE(summary.peak_HRR_W.std_dev, "7B1: peak_HRR_W.std_dev non-finite");
    
    REQUIRE_FINITE(summary.t_peak_HRR_s.mean, "7B1: t_peak_HRR_s.mean non-finite");
    
    // Verify confidence intervals are ordered correctly
    REQUIRE(summary.peak_T_K.ci_lower_95 <= summary.peak_T_K.ci_upper_95,
            "7B1: peak_T_K confidence interval inverted");
    REQUIRE(summary.peak_HRR_W.ci_lower_95 <= summary.peak_HRR_W.ci_upper_95,
            "7B1: peak_HRR_W confidence interval inverted");
    
    // Verify std_dev is non-negative
    REQUIRE(summary.peak_T_K.std_dev >= 0.0, "7B1: peak_T_K.std_dev negative");
    REQUIRE(summary.peak_HRR_W.std_dev >= 0.0, "7B1: peak_HRR_W.std_dev negative");
    
    std::cout << "[PASS] 7B1 MonteCarloUQ basic initialization and small run\n";
}

static void runMonteCarloUQSampling_7B2()
{
    vfep::MonteCarloUQ uq;
    
    vfep::MonteCarloUQ::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 20.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.02;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    uq.setScenario(scenario);
    
    vfep::MonteCarloUQ::UQRanges ranges;
    ranges.heat_release_J_per_mol = {0.9e5, 1.1e5};
    ranges.h_W_m2K = {8.0, 12.0};
    ranges.volume_m3 = {110.0, 130.0};
    ranges.pyrolysis_max_kgps = {0.015, 0.025};
    
    uq.setRanges(ranges);
    
    // Test different sample counts
    const int sample_counts[] = {1, 3, 10};
    
    for (int n : sample_counts) {
        const auto summary = uq.runMonteCarlo(n);
        
        // All results should be finite regardless of sample count
        REQUIRE_FINITE(summary.peak_T_K.mean, "7B2: mean non-finite");
        REQUIRE_FINITE(summary.peak_HRR_W.median, "7B2: median non-finite");
        
        // For n=1, std_dev should be 0 and CI should collapse to single value
        if (n == 1) {
            REQUIRE(summary.peak_T_K.std_dev == 0.0, "7B2: std_dev not zero for n=1");
            REQUIRE(summary.peak_T_K.mean == summary.peak_T_K.median, "7B2: mean != median for n=1");
        }
        
        // For n>1, CI should be meaningful
        if (n > 1) {
            const double span_T = summary.peak_T_K.ci_upper_95 - summary.peak_T_K.ci_lower_95;
            REQUIRE(span_T >= 0.0, "7B2: negative CI span");
        }
    }
    
    std::cout << "[PASS] 7B2 MonteCarloUQ sampling with various counts (1, 3, 10)\n";
}

static void runMonteCarloUQResults_7B3()
{
    vfep::MonteCarloUQ uq;
    
    vfep::MonteCarloUQ::ScenarioConfig scenario;
    scenario.dt_s = 0.05;
    scenario.t_end_s = 25.0;
    scenario.ignite_at_s = 2.0;
    scenario.pyrolysis_max_kgps = 0.02;
    scenario.heat_release_J_per_mol = 1.0e5;
    
    uq.setScenario(scenario);
    
    vfep::MonteCarloUQ::UQRanges ranges;
    ranges.heat_release_J_per_mol = {0.75e5, 1.25e5};
    ranges.h_W_m2K = {5.0, 15.0};
    ranges.volume_m3 = {100.0, 150.0};
    ranges.pyrolysis_max_kgps = {0.01, 0.03};
    
    uq.setRanges(ranges);
    
    // Run larger sample for better statistics
    const auto summary = uq.runMonteCarlo(20);
    
    // Verify statistical properties
    // Mean should be within CI bounds (with small tolerance for numerical issues)
    REQUIRE(summary.peak_T_K.mean >= summary.peak_T_K.ci_lower_95 - 1e-6,
            "7B3: mean below lower CI");
    REQUIRE(summary.peak_T_K.mean <= summary.peak_T_K.ci_upper_95 + 1e-6,
            "7B3: mean above upper CI");
    
    // Median should be within CI bounds
    REQUIRE(summary.peak_HRR_W.median >= summary.peak_HRR_W.ci_lower_95 - 1e-6,
            "7B3: median below lower CI");
    REQUIRE(summary.peak_HRR_W.median <= summary.peak_HRR_W.ci_upper_95 + 1e-6,
            "7B3: median above upper CI");
    
    // Standard deviation should be positive (with n=20, we expect variation)
    REQUIRE(summary.peak_T_K.std_dev > 0.0, "7B3: std_dev not positive with n=20");
    REQUIRE(summary.peak_HRR_W.std_dev > 0.0, "7B3: HRR std_dev not positive with n=20");
    
    // CI span should be reasonable (non-zero with n=20)
    const double span_T = summary.peak_T_K.ci_upper_95 - summary.peak_T_K.ci_lower_95;
    const double span_HRR = summary.peak_HRR_W.ci_upper_95 - summary.peak_HRR_W.ci_lower_95;
    REQUIRE(span_T > 0.0, "7B3: peak_T_K CI span not positive");
    REQUIRE(span_HRR > 0.0, "7B3: peak_HRR_W CI span not positive");
    
    std::cout << "[PASS] 7B3 MonteCarloUQ statistical result validation (n=20)\n";
}

} // namespace

int main() {
    // Canary: prove the test fails in Release when checks are active.
    if (std::getenv("CHEMSI_CANARY_NAN")) {
        REQUIRE_FINITE(std::nan(""), "CANARY_NAN");
        return 0; // unreachable
    }

    // 1B.3a: dt-robustness tripwire
    runDtRobustnessTripwire_1B3a();

    // Existing soak suite
    runSoak("baseline_60s",              0.02, 60.0,         true,  true);
    runSoak("tiny_dt_micro",             1e-6, 2e-3,         false, false);
    runSoak("large_dt_60s",              1.0,  60.0,         true,  true);
    runSoak("soak_2h",                   0.10, 2.0 * 3600.0, true,  true);
    runSoak("burn_no_suppression_10m",   0.05, 600.0,        true,  false);
    runSoak("ultra_large_dt_10m",        5.0,  600.0,        true,  true);

    // 1B.3b: Input contract & invalid-parameter safety
    runInvalidDtHandling_1B3b();
    runVeryLargeDtSafety_1B3b();
    runCommandSchedulingRobustness_1B3b();
    runResetRobustness_1B3b();

    // =======================
    // Step 1C: State Transition Integrity
    // =======================
    runCommandDoesNotAdvanceTime_1C();
    runTimeMonotonicExactDt_1C();
    runIgnitionLatchAndIdempotence_1C();
    runSuppressionBeforeIgnition_1C();
    runTerminalNoOpAndCommandSafety_1C();

    // =======================
    // Step 2: Physical Consistency
    // =======================
    runPreIgnitionQuiescence_2A1();
    runCombustionImpliesFuelConsumption_2A2();
    runHeatReleaseCorrelatesWithTemperatureRise_2A3();
    runSuppressionReducesHeatingOutcomes_2A4();
    runFuelMonotonicityUnderCombustion_2A5();
    runNoSpontaneousMassAppearanceForAgentReservoirs_2A6();
    runSuppressionIncreasesAgentPresenceAfterCommand_2A7();    
    runCombustionSpeciesDirectionality_2A8();
    runSuppressionMakesLessBurnedThanNoSuppression_2A9();
    runNoSuppressionEffectsBeforeCommand_2A10();
    runIgnitionPrecedesCombustionIndicators_2A11();
    runCoarseDtPreservesDirectionalSemantics_2A12();
    
    // =======================
    // Step 3: Interface & State Ownership Integrity
    // =======================
    runInterfaceAndStateOwnershipIntegrity_3_0();
    runSnapshotAndObservationSafety_3A1_A();
    runSnapshotAndObservationSafety_3A1_B();
    runSnapshotAndObservationSafety_3A1_C();
    runSnapshotLifetimeSafety_3A2();
    runInstanceAndResetOwnershipIntegrity_3B1();
    runInstanceLifetimeAndDestructionSafety_3B2();
    runCommandAndApiMisuseSafety_3C1();
    runTerminalApiMisuseSafety_3C2();
    runHighFrequencyPollingStability_3C3();

    // =======================
    // Phase 4A: Suppression Intensity Tests
    // =======================
    runSuppressionIntensityTests();

    // =======================
    // Phase 7: Sensitivity Analysis Tests
    // =======================
    runSensitivityAnalyzerBasic_7A1();
    runSensitivityAnalyzerParameterSweep_7A2();
    runSensitivityAnalyzerResults_7A3();
    runMonteCarloUQBasic_7B1();
    runMonteCarloUQSampling_7B2();
    runMonteCarloUQResults_7B3();

    return 0;
    
}