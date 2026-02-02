#pragma once

#include <vector>

namespace vfep {

class MonteCarloUQ {
public:
    struct ParameterRange {
        double min = 0.0;
        double max = 0.0;
    };

    struct ScenarioGeometry {
        double volume_m3 = 120.0;
        double area_m2 = 180.0;
        double h_W_m2K = 10.0;
    };

    struct ScenarioConfig {
        double dt_s = 0.05;
        double t_end_s = 120.0;
        double ignite_at_s = 2.0;
        double suppress_at_s = 15.0;
        bool enable_suppression = false;
        double ach_1_per_h = -1.0;
        double pyrolysis_max_kgps = 0.03;
        double heat_release_J_per_mol = 1.0e5;
        ScenarioGeometry geometry{};
    };

    struct UQResult {
        double mean = 0.0;
        double median = 0.0;
        double ci_lower_95 = 0.0;
        double ci_upper_95 = 0.0;
        double std_dev = 0.0;
    };

    struct UQSummary {
        UQResult peak_T_K{};
        UQResult peak_HRR_W{};
        UQResult t_peak_HRR_s{};
    };

    struct UQRanges {
        ParameterRange heat_release_J_per_mol{50e3, 200e3};
        ParameterRange h_W_m2K{0.5, 20.0};
        ParameterRange volume_m3{10.0, 500.0};
        ParameterRange pyrolysis_max_kgps{0.01, 1.0};
    };

    MonteCarloUQ();

    void setScenario(const ScenarioConfig& scenario);
    void setRanges(const UQRanges& ranges);

    UQSummary runMonteCarlo(const ScenarioConfig& scenario, int num_samples = 100) const;
    UQSummary runMonteCarlo(int num_samples = 100) const;

private:
    struct SampleMetrics {
        double peak_T_K = 0.0;
        double peak_HRR_W = 0.0;
        double t_peak_HRR_s = 0.0;
    };

    ScenarioConfig scenario_{};
    UQRanges ranges_{};

    SampleMetrics runScenario(const ScenarioConfig& scenario) const;
    UQResult summarize(const std::vector<double>& values) const;
};

} // namespace vfep
