#pragma once
#include <vector>

namespace CloudJ {
class AtmosphericProfile {
    size_t num_layers;
public:
    explicit AtmosphericProfile(size_t layers) : num_layers(layers) {}
    size_t get_num_layers() const noexcept { return num_layers; }
};
}