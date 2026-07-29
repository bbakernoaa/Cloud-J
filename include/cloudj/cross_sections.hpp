#ifndef CLOUDJ_CROSS_SECTIONS_HPP
#define CLOUDJ_CROSS_SECTIONS_HPP

#include <algorithm>

namespace CloudJ {
namespace CrossSections {

/**
 * @brief up-to-three-point interpolation function for X-sections
 * 
 * Supports both standard exact Fortran-matching linear piecewise clamping (default)
 * and branchless GPU-optimized Cubic Hermite Spline (Smoothstep) transition.
 */
inline double interpolate(double t_int, double t1, double x1, double t2, double x2, double t3, double x3, int l123) {
#if defined(CLOUDJ_GPU_MODE)
    // -------------------------------------------------------------------------
    // GPU MODE: 100% branchless, continuously differentiable (C1), SIMD optimal
    // -------------------------------------------------------------------------
    if (l123 <= 1) {
        return x1;
    } else if (l123 == 2) {
        // Normalize and clamp range cleanly
        double t_norm = std::max(0.0, std::min(1.0, (t_int - t1) / (t2 - t1)));
        // Apply Cubic Hermite Spline (Smoothstep: 3*t^2 - 2*t^3)
        double h = t_norm * t_norm * (3.0 - 2.0 * t_norm);
        return x1 + h * (x2 - x1);
    } else {
        // Multi-point spline selection using hardware conditional move select (no branching)
        double t_start = (t_int < t2) ? t1 : t2;
        double t_end   = (t_int < t2) ? t2 : t3;
        double x_start = (t_int < t2) ? x1 : x2;
        double x_end   = (t_int < t2) ? x2 : x3;

        double t_norm = std::max(0.0, std::min(1.0, (t_int - t_start) / (t_end - t_start)));
        double h = t_norm * t_norm * (3.0 - 2.0 * t_norm);
        return x_start + h * (x_end - x_start);
    }
#else
    // -------------------------------------------------------------------------
    // CPU PARITY MODE: Piecewise-linear hard-clamping matching standard Fortran
    // -------------------------------------------------------------------------
    if (l123 <= 1) {
        return x1;
    } else if (l123 == 2) {
        if (t_int <= t1) return x1;
        if (t_int >= t2) return x2;
        return x1 + ((t_int - t1) / (t2 - t1)) * (x2 - x1);
    } else {
        if (t_int <= t1) return x1;
        if (t_int >= t3) return x3;
        if (t_int <= t2) {
            double tfact = (t_int - t1) / (t2 - t1);
            return x1 + tfact * (x2 - x1);
        } else {
            double tfact = (t_int - t2) / (t3 - t2);
            return x2 + tfact * (x3 - x2);
        }
    }
#endif
}

} // namespace CrossSections
} // namespace CloudJ

#endif // CLOUDJ_CROSS_SECTIONS_HPP
