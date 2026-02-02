#pragma once

#include <string>
#include <vector>

namespace vfep {

class SensitivityAnalyzer {
public:
    struct ParameterRange {
        double nominal = 0.0;
        double min = 0.0;
        double max = 0.0;
        int samples = 0;
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

    struct SampleResult {
        double peak_T_K = 0.0;
        double peak_HRR_W = 0.0;
        double t_peak_HRR_s = 0.0;
    };

    struct SensitivityRow {
        std::string parameter_name;
        double parameter_value = 0.0;
        SampleResult metrics{};
    };

    SensitivityAnalyzer();

    void setScenario(const ScenarioConfig& scenario);
    void clearResults();

    void analyzeHeatRelease(const ParameterRange& range);
    void analyzeWallLoss(const ParameterRange& range);
    void analyzeGeometry(const ParameterRange& range);
    void analyzePyrolysis(const ParameterRange& range);

    void exportSensitivityMatrixCSV(const std::string& filename) const;
    const std::vector<SensitivityRow>& results() const;

private:
    ScenarioConfig scenario_{};
    std::vector<SensitivityRow> results_{};

    SampleResult runScenario(const ScenarioConfig& scenario) const;
    std::vector<double> sampleValues(const ParameterRange& range) const;
    static ScenarioGeometry scaleGeometryForVolume(const ScenarioGeometry& base, double volume_m3);
};

} // namespace vfep
