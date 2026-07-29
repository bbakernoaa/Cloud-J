#ifndef CLOUDJ_CROSS_SECTIONS_HPP
#define CLOUDJ_CROSS_SECTIONS_HPP

#include <algorithm>

namespace CloudJ {
namespace CrossSections {

/**
 * @brief up-to-three-point linear interpolation function for X-sections
 * 
 * @param t_int Interpolation temperature (or pressure)
 * @param t1 Temperature for point 1
 * @param x1 Cross-section value for point 1
 * @param t2 Temperature for point 2
 * @param x2 Cross-section value for point 2
 * @param t3 Temperature for point 3
 * @param x3 Cross-section value for point 3
 * @param l123 Number of points to use (1, 2, or 3)
 * @return Interpolated cross-section value
 */
inline double interpolate(double t_int, double t1, double x1, double t2, double x2, double t3, double x3, int l123) {
    if (l123 <= 1) {
        return x1;
    } else if (l123 == 2) {
        double tfact = std::max(0.0, std::min(1.0, (t_int - t1) / (t2 - t1)));
        return x1 + tfact * (x2 - x1);
    } else {
        if (t_int <= t2) {
            double tfact = std::max(0.0, std::min(1.0, (t_int - t1) / (t2 - t1)));
            return x1 + tfact * (x2 - x1);
        } else {
            double tfact = std::max(0.0, std::min(1.0, (t_int - t2) / (t3 - t2)));
            return x2 + tfact * (x3 - x2);
        }
    }
}

} // namespace CrossSections
} // namespace CloudJ

#endif // CLOUDJ_CROSS_SECTIONS_HPP
