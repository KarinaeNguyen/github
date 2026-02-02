#include "UncertaintyQuantification.h"

#include "Simulation.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace vfep {

namespace {
std::vector<double> latinHypercubeSamples(double min_val, double max_val, int samples, std::mt19937& rng) {
    std::vector<double> bins;
    bins.reserve(static_cast<std::size_t>(samples));

    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
    for (int i = 0; i < samples; ++i) {
        const double u = (static_cast<double>(i) + unit_dist(rng)) / static_cast<double>(samples);
        bins.push_back(u);
    }
    std::shuffle(bins.begin(), bins.end(), rng);

    const double span = max_val - min_val;
    for (double& v : bins) {
        v = min_val + span * v;
    }
    return bins;
}

int clampSamples(int samples) {
    return samples < 1 ? 1 : samples;
}

vfep::MonteCarloUQ::ScenarioGeometry scaleGeometryForVolume(
    const vfep::MonteCarloUQ::ScenarioGeometry& base,
    double volume_m3) {
    if (base.volume_m3 <= 0.0 || volume_m3 <= 0.0) {
        return base;
    }

    const double scale = std::cbrt(volume_m3 / base.volume_m3);
    vfep::MonteCarloUQ::ScenarioGeometry scaled = base;
    scaled.volume_m3 = volume_m3;
    scaled.area_m2 = base.area_m2 * (scale * scale);
    return scaled;
}
} // namespace

MonteCarloUQ::MonteCarloUQ() = default;

void MonteCarloUQ::setScenario(const ScenarioConfig& scenario) {
    scenario_ = scenario;
}

void MonteCarloUQ::setRanges(const UQRanges& ranges) {
    ranges_ = ranges;
}

MonteCarloUQ::SampleMetrics MonteCarloUQ::runScenario(const ScenarioConfig& scenario) const {
    vfep::Simulation sim;
    sim.resetToDataCenterRackScenario();

    if (scenario.ach_1_per_h > 0.0) {
        sim.setVentilationACH(scenario.ach_1_per_h);
    }
    sim.setPyrolysisMax(scenario.pyrolysis_max_kgps);
    if (scenario.heat_release_J_per_mol > 0.0) {
        sim.setCombustionHeatRelease(scenario.heat_release_J_per_mol);
    }
    sim.setReactorGeometry(scenario.geometry.volume_m3,
                           scenario.geometry.area_m2,
                           scenario.geometry.h_W_m2K);
    sim.setLiIonEnabled(false);

    SampleMetrics m{};
    double t = 0.0;
    bool ignited = false;
    bool suppressed = false;

    while (t + scenario.dt_s <= scenario.t_end_s + 1e-12) {
        const double t_prev = t;
        const double t_next = t + scenario.dt_s;

        if (!ignited && (t_prev < scenario.ignite_at_s && t_next >= scenario.ignite_at_s)) {
            sim.commandIgniteOrIncreasePyrolysis();
            sim.setPyrolysisRate(scenario.pyrolysis_max_kgps);
            ignited = true;
        }
        if (scenario.enable_suppression && !suppressed &&
            (t_prev < scenario.suppress_at_s && t_next >= scenario.suppress_at_s)) {
            sim.commandStartSuppression();
            sim.setAgentDeliveryRate(1.0);
            sim.setKnockdown(0.55);
            suppressed = true;
        }

        sim.step(scenario.dt_s);
        t = t_next;

        const auto o = sim.observe();
        if (o.T_K > m.peak_T_K) {
            m.peak_T_K = o.T_K;
        }
        if (o.HRR_W > m.peak_HRR_W) {
            m.peak_HRR_W = o.HRR_W;
            m.t_peak_HRR_s = t;
        }
    }

    return m;
}

MonteCarloUQ::UQResult MonteCarloUQ::summarize(const std::vector<double>& values) const {
    UQResult result{};
    if (values.empty()) {
        return result;
    }

    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double variance = 0.0;
    for (double v : values) {
        const double d = v - mean;
        variance += d * d;
    }
    variance /= static_cast<double>(values.size());

    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());

    auto percentile = [&](double p) {
        if (sorted.size() == 1) return sorted.front();
        const double pos = p * (sorted.size() - 1);
        const std::size_t idx = static_cast<std::size_t>(pos);
        const double frac = pos - static_cast<double>(idx);
        if (idx + 1 >= sorted.size()) return sorted.back();
        return sorted[idx] * (1.0 - frac) + sorted[idx + 1] * frac;
    };

    result.mean = mean;
    result.median = percentile(0.5);
    result.ci_lower_95 = percentile(0.025);
    result.ci_upper_95 = percentile(0.975);
    result.std_dev = std::sqrt(variance);
    return result;
}

MonteCarloUQ::UQSummary MonteCarloUQ::runMonteCarlo(const ScenarioConfig& scenario, int num_samples) const {
    const int samples = clampSamples(num_samples);
    std::mt19937 rng(1337u);

    auto heat_samples = latinHypercubeSamples(ranges_.heat_release_J_per_mol.min,
                                              ranges_.heat_release_J_per_mol.max,
                                              samples,
                                              rng);
    auto hW_samples = latinHypercubeSamples(ranges_.h_W_m2K.min,
                                            ranges_.h_W_m2K.max,
                                            samples,
                                            rng);
    auto volume_samples = latinHypercubeSamples(ranges_.volume_m3.min,
                                                ranges_.volume_m3.max,
                                                samples,
                                                rng);
    auto pyro_samples = latinHypercubeSamples(ranges_.pyrolysis_max_kgps.min,
                                              ranges_.pyrolysis_max_kgps.max,
                                              samples,
                                              rng);

    std::vector<double> peak_T;
    std::vector<double> peak_HRR;
    std::vector<double> t_peak_HRR;
    peak_T.reserve(static_cast<std::size_t>(samples));
    peak_HRR.reserve(static_cast<std::size_t>(samples));
    t_peak_HRR.reserve(static_cast<std::size_t>(samples));

    for (int i = 0; i < samples; ++i) {
        ScenarioConfig varied = scenario;
        varied.heat_release_J_per_mol = heat_samples[i];
        varied.geometry.h_W_m2K = hW_samples[i];
        varied.geometry = scaleGeometryForVolume(varied.geometry, volume_samples[i]);
        varied.pyrolysis_max_kgps = pyro_samples[i];

        const auto metrics = runScenario(varied);
        peak_T.push_back(metrics.peak_T_K);
        peak_HRR.push_back(metrics.peak_HRR_W);
        t_peak_HRR.push_back(metrics.t_peak_HRR_s);
    }

    UQSummary summary{};
    summary.peak_T_K = summarize(peak_T);
    summary.peak_HRR_W = summarize(peak_HRR);
    summary.t_peak_HRR_s = summarize(t_peak_HRR);
    return summary;
}

MonteCarloUQ::UQSummary MonteCarloUQ::runMonteCarlo(int num_samples) const {
    return runMonteCarlo(scenario_, num_samples);
}

} // namespace vfep
