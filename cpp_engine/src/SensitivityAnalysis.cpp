#include "SensitivityAnalysis.h"

#include "Simulation.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace vfep {

SensitivityAnalyzer::SensitivityAnalyzer() = default;

void SensitivityAnalyzer::setScenario(const ScenarioConfig& scenario) {
    scenario_ = scenario;
}

void SensitivityAnalyzer::clearResults() {
    results_.clear();
}

std::vector<double> SensitivityAnalyzer::sampleValues(const ParameterRange& range) const {
    std::vector<double> values;
    if (range.samples <= 1 || range.max <= range.min) {
        values.push_back(range.nominal);
        return values;
    }

    values.reserve(static_cast<std::size_t>(range.samples));
    const double span = range.max - range.min;
    const int steps = range.samples - 1;
    for (int i = 0; i < range.samples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(steps);
        values.push_back(range.min + span * t);
    }
    return values;
}

SensitivityAnalyzer::ScenarioGeometry SensitivityAnalyzer::scaleGeometryForVolume(
    const ScenarioGeometry& base,
    double volume_m3) {
    if (base.volume_m3 <= 0.0 || volume_m3 <= 0.0) {
        return base;
    }

    const double scale = std::cbrt(volume_m3 / base.volume_m3);
    ScenarioGeometry scaled = base;
    scaled.volume_m3 = volume_m3;
    scaled.area_m2 = base.area_m2 * (scale * scale);
    return scaled;
}

SensitivityAnalyzer::SampleResult SensitivityAnalyzer::runScenario(const ScenarioConfig& scenario) const {
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

    SampleResult m{};
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

void SensitivityAnalyzer::analyzeHeatRelease(const ParameterRange& range) {
    clearResults();
    for (double value : sampleValues(range)) {
        ScenarioConfig scenario = scenario_;
        scenario.heat_release_J_per_mol = value;
        const auto metrics = runScenario(scenario);
        results_.push_back({"heat_release_J_per_mol", value, metrics});
    }
}

void SensitivityAnalyzer::analyzeWallLoss(const ParameterRange& range) {
    clearResults();
    for (double value : sampleValues(range)) {
        ScenarioConfig scenario = scenario_;
        scenario.geometry.h_W_m2K = value;
        const auto metrics = runScenario(scenario);
        results_.push_back({"h_W_m2K", value, metrics});
    }
}

void SensitivityAnalyzer::analyzeGeometry(const ParameterRange& range) {
    clearResults();
    for (double value : sampleValues(range)) {
        ScenarioConfig scenario = scenario_;
        scenario.geometry = scaleGeometryForVolume(scenario_.geometry, value);
        const auto metrics = runScenario(scenario);
        results_.push_back({"volume_m3", value, metrics});
    }
}

void SensitivityAnalyzer::analyzePyrolysis(const ParameterRange& range) {
    clearResults();
    for (double value : sampleValues(range)) {
        ScenarioConfig scenario = scenario_;
        scenario.pyrolysis_max_kgps = value;
        const auto metrics = runScenario(scenario);
        results_.push_back({"pyrolysis_max_kgps", value, metrics});
    }
}

void SensitivityAnalyzer::exportSensitivityMatrixCSV(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        return;
    }

    out << "parameter,value,peak_T_K,peak_HRR_W,t_peak_HRR_s\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& row : results_) {
        out << row.parameter_name << ','
            << row.parameter_value << ','
            << row.metrics.peak_T_K << ','
            << row.metrics.peak_HRR_W << ','
            << row.metrics.t_peak_HRR_s << '\n';
    }
}

const std::vector<SensitivityAnalyzer::SensitivityRow>& SensitivityAnalyzer::results() const {
    return results_;
}

} // namespace vfep
