#include <cloudj/profile.hpp>
#include <cassert>

void test_profile_bounds() {
    CloudJ::AtmosphericProfile profile(80);
    assert(profile.get_num_layers() == 80);
}

int main() {
    test_profile_bounds();
    return 0;
}
