#ifndef CLOUDJ_PROFILE_HPP
#define CLOUDJ_PROFILE_HPP

#include <vector>

namespace CloudJ {

class AtmosphericProfile {
private:
    size_t num_layers;
    std::vector<double> pressures;     // size num_layers + 1
    std::vector<double> temperatures;  // size num_layers + 1

public:
    explicit AtmosphericProfile(size_t layers) : num_layers(layers) {
        pressures.assign(layers + 1, 0.0);
        temperatures.assign(layers + 1, 0.0);
    }
    
    size_t get_num_layers() const noexcept { return num_layers; }
    
    void set_pressures(const std::vector<double>& p) { pressures = p; }
    void set_temperatures(const std::vector<double>& t) { temperatures = t; }
    
    const std::vector<double>& get_pressures() const noexcept { return pressures; }
    const std::vector<double>& get_temperatures() const noexcept { return temperatures; }
};

} // namespace CloudJ

#endif // CLOUDJ_PROFILE_HPP
