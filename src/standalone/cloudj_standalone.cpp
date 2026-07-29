#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cloudj/cloudj.hpp>

// Convert a double to Fortran's e9.2 scientific formatting cleanly
std::string format_fortran_e9_2(double val) {
    if (std::abs(val) < 1e-99) {
        return " 0.00E+00";
    }
    std::ostringstream ss;
    ss << std::scientific << std::uppercase << std::setprecision(2) << val;
    std::string s = ss.str(); // e.g., "1.23E-05" or "1.23E-005"

    // Parse sign, mantissa, and exponent
    size_t e_pos = s.find('E');
    if (e_pos == std::string::npos) {
        return " 0.00E+00";
    }

    std::string mantissa = s.substr(0, e_pos);
    std::string exp_part = s.substr(e_pos + 1);

    // Parse exponent sign and value
    char exp_sign = '+';
    if (exp_part[0] == '-' || exp_part[0] == '+') {
        exp_sign = exp_part[0];
        exp_part = exp_part.substr(1);
    }

    int exp_val = std::stoi(exp_part);
    std::ostringstream exp_ss;
    exp_ss << exp_sign << std::setw(2) << std::setfill('0') << exp_val;

    // Pad prefix space to fit exactly 9 characters width
    std::string result = mantissa + "E" + exp_ss.str();
    if (result[0] != '-') {
        result = " " + result;
    }
    while (result.length() < 9) {
        result = " " + result;
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::string tables_dir = "tables";
    std::string output_path = "cpp_actual_output.txt";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        }
    }

    // Determine path to profile file dynamically
    std::string profile_filename = tables_dir + "/atmos_PTClds.dat";
    std::ifstream infile(profile_filename);
    if (!infile.is_open()) {
        // Fallback to look under build/bin/tables if executed inside CMake bin
        profile_filename = "./tables/atmos_PTClds.dat";
        infile.open(profile_filename);
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open profile file tables/atmos_PTClds.dat\n";
            return 1;
        }
    }

    std::string line;
    // Skip first 2 lines
    std::getline(infile, line);
    std::getline(infile, line);

    // Read surface pressure (line 3)
    double psurf = 1013.25;
    if (std::getline(infile, line)) {
        std::stringstream ss(line);
        ss >> psurf;
    }

    // Skip lines 4, 5, 6, 7
    for (int i = 0; i < 4; ++i) {
        std::getline(infile, line);
    }

    // Read 58 levels of data
    std::vector<double> etaa(58);
    std::vector<double> etab(58);
    std::vector<double> temp(58);
    std::vector<double> rh(58);

    for (int i = 0; i < 58; ++i) {
        if (std::getline(infile, line)) {
            std::stringstream ss(line);
            int idx;
            ss >> idx >> etaa[i] >> etab[i] >> temp[i] >> rh[i];
        }
    }
    infile.close();

    // Calculate boundary/edges pressures and temperatures
    std::vector<double> pressures(58);
    std::vector<double> temperatures(58);
    for (int i = 0; i < 58; ++i) {
        pressures[i] = etaa[i] + etab[i] * psurf;
        temperatures[i] = temp[i];
    }

    CloudJ::Engine engine;
    
    // We instantiate AtmosphericProfile with 57 layers (which has 58 boundary nodes/levels)
    CloudJ::AtmosphericProfile profile(57);
    profile.set_pressures(pressures);
    profile.set_temperatures(temperatures);

    std::cout << ">>>begin Cloud-J v8.0 C++ Standalone\n";

    // Replicate Fortran Standalone SZA scan for index loop 1 (SZA = 0)
    // SZA values: SZAscan = [0, 30, 60]
    std::vector<double> SZAscan = {0.0, 30.0, 60.0};
    
    for (size_t i = 0; i < 3; ++i) {
        double sza = SZAscan[i];
        CloudJ::OutputRates rates = engine.calculate_photolysis_rates(profile, sza);

        std::cout << " Fast-J ----J-values----\n";
        std::cout << "L=  ";
        for (const auto& title : engine.get_spec_data().titlejx) {
            std::cout << std::setw(6) << title << "   ";
        }
        std::cout << "\n";

        // Print levels 57 down to 1 (0-based 56 down to 0)
        for (int l = 57; l >= 1; --l) {
            std::cout << std::setw(3) << l;
            for (int j = 0; j < engine.get_spec_data().njx; ++j) {
                std::cout << format_fortran_e9_2(rates.j_values[l - 1][j]);
            }
            std::cout << "\n";
        }
    }

    std::cout << "C++ standalone calculation completed successfully!\n";
    return 0;
}
