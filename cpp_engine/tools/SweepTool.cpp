#include "SensitivityAnalysis.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>

namespace {
std::string toLower(std::string v) {
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return v;
}

void printUsage() {
    std::cout << "SweepTool usage:\n"
              << "  SweepTool --param <heat_release|h_w|volume|pyrolysis> [--min v] [--max v] [--samples n] [--out file]\n";
}
} // namespace

int main(int argc, char** argv) {
    std::string param;
    double min_val = 0.0;
    double max_val = 0.0;
    int samples = 5;
    std::string out = "sensitivity.csv";
    bool min_set = false;
    bool max_set = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--param" && i + 1 < argc) {
            param = toLower(argv[++i]);
        } else if (arg == "--min" && i + 1 < argc) {
            min_val = std::stod(argv[++i]);
            min_set = true;
        } else if (arg == "--max" && i + 1 < argc) {
            max_val = std::stod(argv[++i]);
            max_set = true;
        } else if (arg == "--samples" && i + 1 < argc) {
            samples = std::stoi(argv[++i]);
        } else if (arg == "--out" && i + 1 < argc) {
            out = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        } else {
            std::cout << "Unknown argument: " << arg << "\n";
            printUsage();
            return 1;
        }
    }

    if (param.empty()) {
        printUsage();
        return 1;
    }

    vfep::SensitivityAnalyzer analyzer;
    vfep::SensitivityAnalyzer::ScenarioConfig scenario;
    analyzer.setScenario(scenario);

    vfep::SensitivityAnalyzer::ParameterRange range;
    range.samples = samples;

    if (param == "heat_release" || param == "heat_release_j_per_mol") {
        range.nominal = scenario.heat_release_J_per_mol;
    } else if (param == "h_w" || param == "h_w_m2k" || param == "wall_loss") {
        range.nominal = scenario.geometry.h_W_m2K;
    } else if (param == "volume" || param == "volume_m3") {
        range.nominal = scenario.geometry.volume_m3;
    } else if (param == "pyrolysis" || param == "pyrolysis_max" || param == "pyrolysis_max_kgps") {
        range.nominal = scenario.pyrolysis_max_kgps;
    } else {
        std::cout << "Unsupported parameter: " << param << "\n";
        printUsage();
        return 1;
    }

    if (!min_set) {
        min_val = range.nominal * 0.75;
    }
    if (!max_set) {
        max_val = range.nominal * 1.25;
    }

    range.min = min_val;
    range.max = max_val;

    if (param == "heat_release" || param == "heat_release_j_per_mol") {
        analyzer.analyzeHeatRelease(range);
    } else if (param == "h_w" || param == "h_w_m2k" || param == "wall_loss") {
        analyzer.analyzeWallLoss(range);
    } else if (param == "volume" || param == "volume_m3") {
        analyzer.analyzeGeometry(range);
    } else {
        analyzer.analyzePyrolysis(range);
    }

    analyzer.exportSensitivityMatrixCSV(out);
    std::cout << "Wrote sensitivity sweep to: " << out << "\n";
    return 0;
}
