#ifndef CLOUDJ_RADIATIVE_SOLVER_HPP
#define CLOUDJ_RADIATIVE_SOLVER_HPP

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <experimental/mdspan.hpp>

// Portable compiler loop-unrolling hint macros for standard compilers (GCC, Clang, Intel oneAPI)
#if defined(__clang__)
  #define CLOUDJ_UNROLL_4 _Pragma("clang loop unroll_count(4)")
#elif defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER)
  #define CLOUDJ_UNROLL_4 _Pragma("unroll(4)")
#elif defined(__GNUC__)
  #define CLOUDJ_UNROLL_4 _Pragma("GCC unroll 4")
#else
  #define CLOUDJ_UNROLL_4
#endif

namespace CloudJ {
namespace RadiativeSolver {

constexpr int M_ = 4;
constexpr int M2_ = 8;

struct alignas(64) Workspace {
    // 2D buffers flat storage
    std::vector<double> a_data;  // size: M_ * nd
    std::vector<double> c_data;  // size: M_ * nd
    std::vector<double> h_data;  // size: M_ * nd
    std::vector<double> rr_data; // size: M_ * nd

    // 3D buffers flat storage
    std::vector<double> b_data;  // size: M_ * M_ * nd
    std::vector<double> aa_data; // size: M_ * M_ * nd
    std::vector<double> cc_data; // size: M_ * M_ * nd
    std::vector<double> dd_data; // size: M_ * M_ * nd

    // Resizes all vectors once to the required layer depth nd
    void resize(size_t nd) {
        size_t size_2d = M_ * nd;
        size_t size_3d = M_ * M_ * nd;

        a_data.assign(size_2d, 0.0);
        c_data.assign(size_2d, 0.0);
        h_data.assign(size_2d, 0.0);
        rr_data.assign(size_2d, 0.0);

        b_data.assign(size_3d, 0.0);
        aa_data.assign(size_3d, 0.0);
        cc_data.assign(size_3d, 0.0);
        dd_data.assign(size_3d, 0.0);
    }
};

// standard 8-stream Gauss points & weights from cldj_cmn_mod.F90
constexpr std::array<double, M_> EMU = {
    0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703
};
constexpr std::array<double, M_> WT = {
    0.17392742256873, 0.32607257743127, 0.32607257743127, 0.17392742256873
};

using mdspan_2d = std::experimental::mdspan<const double, std::experimental::dextents<size_t, 2>, std::experimental::layout_left>;
using mdspan_1d = std::experimental::mdspan<const double, std::experimental::dextents<size_t, 1>, std::experimental::layout_left>;

using mdspan_3d_mut = std::experimental::mdspan<double, std::experimental::dextents<size_t, 3>, std::experimental::layout_left>;
using mdspan_2d_mut = std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>, std::experimental::layout_left>;

/**
 * @brief Generates coefficient matrices for the block tri-diagonal system.
 * Matches GEN_ID in cldj_fjx_sub_mod.F90.
 */
inline void GEN_ID(
    mdspan_2d pomega,  // (M2_, N_)
    mdspan_1d fz,      // (N_)
    mdspan_1d ztau,    // (N_)
    double zflux,
    const std::array<double, 5>& rfl,
    const double pm[M_][M2_],
    const double pm0[M2_],
    mdspan_3d_mut b,   // (M_, M_, N_)
    mdspan_3d_mut aa,  // (M_, M_, N_)
    mdspan_3d_mut cc,  // (M_, M_, N_)
    mdspan_2d_mut a,   // (M_, N_)
    mdspan_2d_mut h,   // (M_, N_)
    mdspan_2d_mut c,   // (M_, N_)
    int nd
) {
    // Local 4x4 matrix helpers
    double s[M_][M_] = {0};
    double t[M_][M_] = {0};
    double u[M_][M_] = {0};
    double v[M_][M_] = {0};
    double w[M_][M_] = {0};

    // Initialize outputs
    for (int l = 0; l < nd; ++l) {
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
            a(i, l) = 0.0;
            c(i, l) = 0.0;
            h(i, l) = 0.0;
            CLOUDJ_UNROLL_4
            for (int j = 0; j < M_; ++j) {
                b(i, j, l) = 0.0;
                aa(i, j, l) = 0.0;
                cc(i, j, l) = 0.0;
            }
        }
    }

    // Upper boundary: 2nd-order terms
    int l1 = 0; // 0-based Fortran L=1
    int l2 = 1; // 0-based Fortran L=2

    for (int i = 0; i < M_; ++i) {
        double sum0 = pomega(0, l1) * pm[i][0] * pm0[0] + pomega(2, l1) * pm[i][2] * pm0[2] +
                      pomega(4, l1) * pm[i][4] * pm0[4] + pomega(6, l1) * pm[i][6] * pm0[6];
        double sum2 = pomega(0, l2) * pm[i][0] * pm0[0] + pomega(2, l2) * pm[i][2] * pm0[2] +
                      pomega(4, l2) * pm[i][4] * pm0[4] + pomega(6, l2) * pm[i][6] * pm0[6];
        double sum1 = pomega(1, l1) * pm[i][1] * pm0[1] + pomega(3, l1) * pm[i][3] * pm0[3] +
                      pomega(5, l1) * pm[i][5] * pm0[5] + pomega(7, l1) * pm[i][7] * pm0[7];
        double sum3 = pomega(1, l2) * pm[i][1] * pm0[1] + pomega(3, l2) * pm[i][3] * pm0[3] +
                      pomega(5, l2) * pm[i][5] * pm0[5] + pomega(7, l2) * pm[i][7] * pm0[7];
        h(i, l1) = 0.5 * (sum0 * fz(l1) + sum2 * fz(l2));
        a(i, l1) = 0.5 * (sum1 * fz(l1) + sum3 * fz(l2));
    }

    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum0 = pomega(0, l1) * pm[i][0] * pm[j][0] + pomega(2, l1) * pm[i][2] * pm[j][2] +
                          pomega(4, l1) * pm[i][4] * pm[j][4] + pomega(6, l1) * pm[i][6] * pm[j][6];
            double sum2 = pomega(0, l2) * pm[i][0] * pm[j][0] + pomega(2, l2) * pm[i][2] * pm[j][2] +
                          pomega(4, l2) * pm[i][4] * pm[j][4] + pomega(6, l2) * pm[i][6] * pm[j][6];
            double sum1 = pomega(1, l1) * pm[i][1] * pm[j][1] + pomega(3, l1) * pm[i][3] * pm[j][3] +
                          pomega(5, l1) * pm[i][5] * pm[j][5] + pomega(7, l1) * pm[i][7] * pm[j][7];
            double sum3 = pomega(1, l2) * pm[i][1] * pm[j][1] + pomega(3, l2) * pm[i][3] * pm[j][3] +
                          pomega(5, l2) * pm[i][5] * pm[j][5] + pomega(7, l2) * pm[i][7] * pm[j][7];
            
            s[i][j] = -sum2 * WT[j];
            s[j][i] = -sum2 * WT[i];
            t[i][j] = -sum1 * WT[j];
            t[j][i] = -sum1 * WT[i];
            v[i][j] = -sum3 * WT[j];
            v[j][i] = -sum3 * WT[i];
            b(i, j, l1) = -0.5 * (sum0 + sum2) * WT[j];
            b(j, i, l1) = -0.5 * (sum0 + sum2) * WT[i];
        }
    }

    for (int i = 0; i < M_; ++i) {
        s[i][i] += 1.0;
        t[i][i] += 1.0;
        v[i][i] += 1.0;
        b(i, i, l1) += 1.0;

        c(i, l1) = s[i][0] * a(0, l1) / EMU[0] + s[i][1] * a(1, l1) / EMU[1] +
                   s[i][2] * a(2, l1) / EMU[2] + s[i][3] * a(3, l1) / EMU[3];
    }

    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            w[j][i] = s[j][0] * t[0][i] / EMU[0] + s[j][1] * t[1][i] / EMU[1] +
                      s[j][2] * t[2][i] / EMU[2] + s[j][3] * t[3][i] / EMU[3];
            u[j][i] = s[j][0] * v[0][i] / EMU[0] + s[j][1] * v[1][i] / EMU[1] +
                      s[j][2] * v[2][i] / EMU[2] + s[j][3] * v[3][i] / EMU[3];
        }
    }

    double deltau = ztau(l2) - ztau(l1);
    double d2 = 0.25 * deltau;
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            b(i, j, l1) += d2 * w[i][j];
            cc(i, j, l1) = d2 * u[i][j];
        }
        h(i, l1) += 2.0 * d2 * c(i, l1);
        a(i, l1) = 0.0;
    }

    for (int i = 0; i < M_; ++i) {
        double d1 = EMU[i] / deltau;
        b(i, i, l1) += d1;
        cc(i, i, l1) -= d1;
    }

    // Intermediate points: can be even or odd, A & C diagonal
    for (int ll = 2; ll <= nd - 2; ll += 2) {
        deltau = ztau(ll + 1) - ztau(ll - 1);
        for (int i = 0; i < M_; ++i) {
            a(i, ll) = EMU[i] / deltau;
            c(i, ll) = -a(i, ll);
            h(i, ll) = fz(ll) * (pomega(1, ll) * pm[i][1] * pm0[1] + pomega(3, ll) * pm[i][3] * pm0[3] +
                                 pomega(5, ll) * pm[i][5] * pm0[5] + pomega(7, ll) * pm[i][7] * pm0[7]);
        }
        for (int i = 0; i < M_; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum0 = pomega(1, ll) * pm[i][1] * pm[j][1] + pomega(3, ll) * pm[i][3] * pm[j][3] +
                              pomega(5, ll) * pm[i][5] * pm[j][5] + pomega(7, ll) * pm[i][7] * pm[j][7];
                b(i, j, ll) = -sum0 * WT[j];
                b(j, i, ll) = -sum0 * WT[i];
            }
        }
        for (int i = 0; i < M_; ++i) {
            b(i, i, ll) += 1.0;
        }
    }

    for (int ll = 3; ll <= nd - 3; ll += 2) {
        deltau = ztau(ll + 1) - ztau(ll - 1);
        for (int i = 0; i < M_; ++i) {
            a(i, ll) = EMU[i] / deltau;
            c(i, ll) = -a(i, ll);
            h(i, ll) = fz(ll) * (pomega(0, ll) * pm[i][0] * pm0[0] + pomega(2, ll) * pm[i][2] * pm0[2] +
                                 pomega(4, ll) * pm[i][4] * pm0[4] + pomega(6, ll) * pm[i][6] * pm0[6]);
        }
        for (int i = 0; i < M_; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum0 = pomega(0, ll) * pm[i][0] * pm[j][0] + pomega(2, ll) * pm[i][2] * pm[j][2] +
                              pomega(4, ll) * pm[i][4] * pm[j][4] + pomega(6, ll) * pm[i][6] * pm[j][6];
                b(i, j, ll) = -sum0 * WT[j];
                b(j, i, ll) = -sum0 * WT[i];
            }
        }
        for (int i = 0; i < M_; ++i) {
            b(i, i, ll) += 1.0;
        }
    }

    // Lower boundary: 2nd-order terms
    int l_last = nd - 1; // 0-based Fortran L=ND
    int l_prev = nd - 2; // 0-based Fortran L=ND-1

    for (int i = 0; i < M_; ++i) {
        double sum0 = pomega(0, l_last) * pm[i][0] * pm0[0] + pomega(2, l_last) * pm[i][2] * pm0[2] +
                      pomega(4, l_last) * pm[i][4] * pm0[4] + pomega(6, l_last) * pm[i][6] * pm0[6];
        double sum2 = pomega(0, l_prev) * pm[i][0] * pm0[0] + pomega(2, l_prev) * pm[i][2] * pm0[2] +
                      pomega(4, l_prev) * pm[i][4] * pm0[4] + pomega(6, l_prev) * pm[i][6] * pm0[6];
        double sum1 = pomega(1, l_last) * pm[i][1] * pm0[1] + pomega(3, l_last) * pm[i][3] * pm0[3] +
                      pomega(5, l_last) * pm[i][5] * pm0[5] + pomega(7, l_last) * pm[i][7] * pm0[7];
        double sum3 = pomega(1, l_prev) * pm[i][1] * pm0[1] + pomega(3, l_prev) * pm[i][3] * pm0[3] +
                      pomega(5, l_prev) * pm[i][5] * pm0[5] + pomega(7, l_prev) * pm[i][7] * pm0[7];
        h(i, l_last) = 0.5 * (sum0 * fz(l_last) + sum2 * fz(l_prev));
        a(i, l_last) = 0.5 * (sum1 * fz(l_last) + sum3 * fz(l_prev));
    }

    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum0 = pomega(0, l_last) * pm[i][0] * pm[j][0] + pomega(2, l_last) * pm[i][2] * pm[j][2] +
                          pomega(4, l_last) * pm[i][4] * pm[j][4] + pomega(6, l_last) * pm[i][6] * pm[j][6];
            double sum2 = pomega(0, l_prev) * pm[i][0] * pm[j][0] + pomega(2, l_prev) * pm[i][2] * pm[j][2] +
                          pomega(4, l_prev) * pm[i][4] * pm[j][4] + pomega(6, l_prev) * pm[i][6] * pm[j][6];
            double sum1 = pomega(1, l_last) * pm[i][1] * pm[j][1] + pomega(3, l_last) * pm[i][3] * pm[j][3] +
                          pomega(5, l_last) * pm[i][5] * pm[j][5] + pomega(7, l_last) * pm[i][7] * pm[j][7];
            double sum3 = pomega(1, l_prev) * pm[i][1] * pm[j][1] + pomega(3, l_prev) * pm[i][3] * pm[j][3] +
                          pomega(5, l_prev) * pm[i][5] * pm[j][5] + pomega(7, l_prev) * pm[i][7] * pm[j][7];
            s[i][j] = -sum2 * WT[j];
            s[j][i] = -sum2 * WT[i];
            t[i][j] = -sum1 * WT[j];
            t[j][i] = -sum1 * WT[i];
            v[i][j] = -sum3 * WT[j];
            v[j][i] = -sum3 * WT[i];
            b(i, j, l_last) = -0.5 * (sum0 + sum2) * WT[j];
            b(j, i, l_last) = -0.5 * (sum0 + sum2) * WT[i];
        }
    }

    for (int i = 0; i < M_; ++i) {
        s[i][i] += 1.0;
        t[i][i] += 1.0;
        v[i][i] += 1.0;
        b(i, i, l_last) += 1.0;

        c(i, l_last) = s[i][0] * a(0, l_last) / EMU[0] + s[i][1] * a(1, l_last) / EMU[1] +
                       s[i][2] * a(2, l_last) / EMU[2] + s[i][3] * a(3, l_last) / EMU[3];
    }

    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            w[j][i] = s[j][0] * t[0][i] / EMU[0] + s[j][1] * t[1][i] / EMU[1] +
                      s[j][2] * t[2][i] / EMU[2] + s[j][3] * t[3][i] / EMU[3];
            u[j][i] = s[j][0] * v[0][i] / EMU[0] + s[j][1] * v[1][i] / EMU[1] +
                      s[j][2] * v[2][i] / EMU[2] + s[j][3] * v[3][i] / EMU[3];
        }
    }

    deltau = ztau(l_last) - ztau(l_prev);
    d2 = 0.25 * deltau;
    double sumrfl = 0.0;
    for (int j = 0; j < M_; ++j) {
        sumrfl += rfl[j] * EMU[j] * WT[j];
    }
    double surfac = 4.0 / (1.0 + 2.0 * sumrfl);

    for (int i = 0; i < M_; ++i) {
        double d1 = EMU[i] / deltau;
        double sum0 = d1 + d2 * (w[i][0] + w[i][1] + w[i][2] + w[i][3]);
        for (int j = 0; j < M_; ++j) {
            aa(i, j, l_last) = -d2 * u[i][j];
            b(i, j, l_last) = b(i, j, l_last) + d2 * w[i][j] - sum0 * surfac * rfl[j] * EMU[j] * WT[j];
        }
        h(i, l_last) = h(i, l_last) - 2.0 * d2 * c(i, l_last) + sum0 * surfac * 0.25 * rfl[4] * zflux;
    }

    for (int i = 0; i < M_; ++i) {
        double d1 = EMU[i] / deltau;
        aa(i, i, l_last) += d1;
        b(i, i, l_last) += d1;
        c(i, l_last) = 0.0;
    }
}

// 4x4 LU solver helper matching manual BLKSLV algorithm (Highly Optimized Reciprocal Form)
inline void solve_lu_4x4(double E[M_][M_]) {
    double inv_E00 = 1.0 / E[0][0];
    E[1][0] *= inv_E00;
    E[1][1] = E[1][1] - E[1][0] * E[0][1];
    E[1][2] = E[1][2] - E[1][0] * E[0][2];
    E[1][3] = E[1][3] - E[1][0] * E[0][3];
    
    E[2][0] *= inv_E00;
    double inv_E11 = 1.0 / E[1][1];
    E[2][1] = (E[2][1] - E[2][0] * E[0][1]) * inv_E11;
    E[2][2] = E[2][2] - E[2][0] * E[0][2] - E[2][1] * E[1][2];
    E[2][3] = E[2][3] - E[2][0] * E[0][3] - E[2][1] * E[1][3];
    
    E[3][0] *= inv_E00;
    E[3][1] = (E[3][1] - E[3][0] * E[0][1]) * inv_E11;
    double inv_E22 = 1.0 / E[2][2];
    E[3][2] = (E[3][2] - E[3][0] * E[0][2] - E[3][1] * E[1][2]) * inv_E22;
    E[3][3] = E[3][3] - E[3][0] * E[0][3] - E[3][1] * E[1][3] - E[3][2] * E[2][3];

    // Invert L
    E[3][2] = -E[3][2];
    E[3][1] = -E[3][1] - E[3][2] * E[2][1];
    E[3][0] = -E[3][0] - E[3][1] * E[1][0] - E[3][2] * E[2][0];
    E[2][1] = -E[2][1];
    E[2][0] = -E[2][0] - E[2][1] * E[1][0];
    E[1][0] = -E[1][0];

    // Invert U (Using pre-calculated reciprocals for division-free speedups)
    E[3][3] = 1.0 / E[3][3];
    E[2][3] = -E[2][3] * E[3][3] * inv_E22;
    E[2][2] = inv_E22;
    E[1][3] = -(E[1][2] * E[2][3] + E[1][3] * E[3][3]) * inv_E11;
    E[1][2] = -E[1][2] * E[2][2] * inv_E11;
    E[1][1] = inv_E11;
    E[0][3] = -(E[0][1] * E[1][3] + E[0][2] * E[2][3] + E[0][3] * E[3][3]) * inv_E00;
    E[0][2] = -(E[0][1] * E[1][2] + E[0][2] * E[2][2]) * inv_E00;
    E[0][1] = -E[0][1] * E[1][1] * inv_E00;
    E[0][0] = inv_E00;

    // Multiply U-inverse * L-inverse, storing result in E
    double temp[M_][M_];
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            temp[i][j] = E[i][j];
        }
    }

    E[0][0] = temp[0][0] + temp[0][1] * temp[1][0] + temp[0][2] * temp[2][0] + temp[0][3] * temp[3][0];
    E[0][1] = temp[0][1] + temp[0][2] * temp[2][1] + temp[0][3] * temp[3][1];
    E[0][2] = temp[0][2] + temp[0][3] * temp[3][2];
    E[1][0] = temp[1][1] * temp[1][0] + temp[1][2] * temp[2][0] + temp[1][3] * temp[3][0];
    E[1][1] = temp[1][1] + temp[1][2] * temp[2][1] + temp[1][3] * temp[3][1];
    E[1][2] = temp[1][2] + temp[1][3] * temp[3][2];
    E[2][0] = temp[2][2] * temp[2][0] + temp[2][3] * temp[3][0];
    E[2][1] = temp[2][2] * temp[2][1] + temp[2][3] * temp[3][1];
    E[2][2] = temp[2][2] + temp[2][3] * temp[3][2];
    E[3][0] = temp[3][3] * temp[3][0];
    E[3][1] = temp[3][3] * temp[3][1];
    E[3][2] = temp[3][3] * temp[3][2];
}

#if defined(CLOUDJ_USE_PCR)
// High-performance 4x4 matrix multiplication helper
inline void mat_mult_4x4(const double X[M_][M_], const double Y[M_][M_], double Z[M_][M_]) {
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
        CLOUDJ_UNROLL_4
        for (int j = 0; j < M_; ++j) {
            Z[i][j] = X[i][0] * Y[0][j] + X[i][1] * Y[1][j] + X[i][2] * Y[2][j] + X[i][3] * Y[3][j];
        }
    }
}

// High-performance 4x4 matrix by 4-vector multiplication helper
inline void mat_vec_mult_4(const double X[M_][M_], const double V[M_], double R[M_]) {
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
        R[i] = X[i][0] * V[0] + X[i][1] * V[1] + X[i][2] * V[2] + X[i][3] * V[3];
    }
}

// Inverts a 4x4 matrix in-place using our optimized division-free LU solver
inline void invert_matrix_4x4_helper(const double B[M_][M_], double invB[M_][M_]) {
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
        CLOUDJ_UNROLL_4
        for (int j = 0; j < M_; ++j) {
            invB[i][j] = B[i][j];
        }
    }
    solve_lu_4x4(invB);
}

inline void solve_pcr(
    mdspan_2d_mut fj,
    mdspan_2d pomega,
    mdspan_1d fz,
    mdspan_1d ztau,
    double fsbot,
    const std::array<double, 5>& rfl,
    const double pm[M_][M2_],
    const double pm0[M2_],
    double& fjtop,
    double& fjbot,
    std::array<double, 5>& fibot,
    int nd,
    int k_idx,
    Workspace& ws
) {
    // Create temporary block-tridiagonal views to assemble the global system coefficients
    mdspan_2d_mut a(ws.a_data.data(), M_, nd);
    mdspan_2d_mut c(ws.c_data.data(), M_, nd);
    mdspan_2d_mut h(ws.h_data.data(), M_, nd);
    mdspan_2d_mut rr(ws.rr_data.data(), M_, nd);

    mdspan_3d_mut b(ws.b_data.data(), M_, M_, nd);
    mdspan_3d_mut aa(ws.aa_data.data(), M_, M_, nd);
    mdspan_3d_mut cc(ws.cc_data.data(), M_, M_, nd);
    mdspan_3d_mut dd(ws.dd_data.data(), M_, M_, nd);

    // Generate block tri-diagonal system coefficients (a, b, cc, c)
    GEN_ID(pomega, fz, ztau, fsbot, rfl, pm, pm0, b, aa, cc, a, h, c, nd);

    // Setup persistent Parallel Cyclic Reduction coefficient scratchpad vectors
    // To ensure zero allocation at runtime, we can utilize the dd_data and b_data workspace pools 
    // or allocate small temporary buffers since this is CPU-forced testing.
    std::vector<double> A_pcr(M_ * M_ * nd, 0.0);
    std::vector<double> B_pcr(M_ * M_ * nd, 0.0);
    std::vector<double> C_pcr(M_ * M_ * nd, 0.0);
    std::vector<double> D_pcr(M_ * nd, 0.0);

    mdspan_3d_mut A_v(A_pcr.data(), M_, M_, nd);
    mdspan_3d_mut B_v(B_pcr.data(), M_, M_, nd);
    mdspan_3d_mut C_v(C_pcr.data(), M_, M_, nd);
    mdspan_2d_mut D_v(D_pcr.data(), M_, nd);

    // Initialize PCR matrices
    for (int l = 0; l < nd; ++l) {
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
            D_v(i, l) = c(i, l);
            // a is diagonal of size 4 representing lower diagonal
            A_v(i, i, l) = a(i, l);
            CLOUDJ_UNROLL_4
            for (int j = 0; j < M_; ++j) {
                B_v(i, j, l) = b(i, j, l);
                C_v(i, j, l) = cc(i, j, l);
            }
        }
    }

    // Boundary conditions adjustments for Thomas equivalence
    // Thomas boundary elements are already generated in GEN_ID. We adjust lower boundary directly.
    int l_last = nd - 1;
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            A_v(i, j, l_last) = aa(i, j, l_last);
        }
    }

    // Execute O(log N) stride-reduction stages
    int stages = std::ceil(std::log2(nd));
    for (int step = 0; step < stages; ++step) {
        int stride = 1 << step;

        std::vector<double> A_next(M_ * M_ * nd, 0.0);
        std::vector<double> B_next(M_ * M_ * nd, 0.0);
        std::vector<double> C_next(M_ * M_ * nd, 0.0);
        std::vector<double> D_next(M_ * nd, 0.0);

        mdspan_3d_mut An(A_next.data(), M_, M_, nd);
        mdspan_3d_mut Bn(B_next.data(), M_, M_, nd);
        mdspan_3d_mut Cn(C_next.data(), M_, M_, nd);
        mdspan_2d_mut Dn(D_next.data(), M_, nd);

        for (int l = 0; l < nd; ++l) {
            double alpha[M_][M_] = {0};
            double beta[M_][M_]  = {0};

            // 1. Calculate alpha = - A_l * B_{l-stride}^-1
            if (l - stride >= 0) {
                double Bl_left[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Bl_left[i][j] = B_v(i, j, l - stride);
                    }
                }
                double invBl_left[M_][M_];
                invert_matrix_4x4_helper(Bl_left, invBl_left);

                double Al[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Al[i][j] = A_v(i, j, l);
                    }
                }
                double temp[M_][M_];
                mat_mult_4x4(Al, invBl_left, temp);
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        alpha[i][j] = -temp[i][j];
                    }
                }
            }

            // 2. Calculate beta = - C_l * B_{l+stride}^-1
            if (l + stride < nd) {
                double Bl_right[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Bl_right[i][j] = B_v(i, j, l + stride);
                    }
                }
                double invBl_right[M_][M_];
                invert_matrix_4x4_helper(Bl_right, invBl_right);

                double Cl[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Cl[i][j] = C_v(i, j, l);
                    }
                }
                double temp[M_][M_];
                mat_mult_4x4(Cl, invBl_right, temp);
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        beta[i][j] = -temp[i][j];
                    }
                }
            }

            // 3. Compute new coefficients (An, Bn, Cn, Dn)
            // Bn = B_l + alpha * C_{l-stride} + beta * A_{l+stride}
            double term1[M_][M_] = {0};
            double term2[M_][M_] = {0};
            if (l - stride >= 0) {
                double Cl_left[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Cl_left[i][j] = C_v(i, j, l - stride);
                    }
                }
                mat_mult_4x4(alpha, Cl_left, term1);
            }
            if (l + stride < nd) {
                double Al_right[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Al_right[i][j] = A_v(i, j, l + stride);
                    }
                }
                mat_mult_4x4(beta, Al_right, term2);
            }

            CLOUDJ_UNROLL_4
            for (int i = 0; i < M_; ++i) {
                CLOUDJ_UNROLL_4
                for (int j = 0; j < M_; ++j) {
                    Bn(i, j, l) = B_v(i, j, l) + term1[i][j] + term2[i][j];
                }
            }

            // An = alpha * A_{l-stride}
            if (l - stride >= 0) {
                double Al_left[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Al_left[i][j] = A_v(i, j, l - stride);
                    }
                }
                double temp[M_][M_];
                mat_mult_4x4(alpha, Al_left, temp);
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        An(i, j, l) = temp[i][j];
                    }
                }
            }

            // Cn = beta * C_{l+stride}
            if (l + stride < nd) {
                double Cl_right[M_][M_];
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Cl_right[i][j] = C_v(i, j, l + stride);
                    }
                }
                double temp[M_][M_];
                mat_mult_4x4(beta, Cl_right, temp);
                CLOUDJ_UNROLL_4
                for (int i = 0; i < M_; ++i) {
                    CLOUDJ_UNROLL_4
                    for (int j = 0; j < M_; ++j) {
                        Cn(i, j, l) = temp[i][j];
                    }
                }
            }

            // Dn = D_l + alpha * D_{l-stride} + beta * D_{l+stride}
            double termD1[M_] = {0};
            double termD2[M_] = {0};
            if (l - stride >= 0) {
                double Dl_left[M_];
                for (int i = 0; i < M_; ++i) Dl_left[i] = D_v(i, l - stride);
                mat_vec_mult_4(alpha, Dl_left, termD1);
            }
            if (l + stride < nd) {
                double Dl_right[M_];
                for (int i = 0; i < M_; ++i) Dl_right[i] = D_v(i, l + stride);
                mat_vec_mult_4(beta, Dl_right, termD2);
            }
            CLOUDJ_UNROLL_4
            for (int i = 0; i < M_; ++i) {
                Dn(i, l) = D_v(i, l) + termD1[i] + termD2[i];
            }
        }

        // Copy new coefficients back to PCR views for next stage
        for (int l = 0; l < nd; ++l) {
            CLOUDJ_UNROLL_4
            for (int i = 0; i < M_; ++i) {
                D_v(i, l) = Dn(i, l);
                CLOUDJ_UNROLL_4
                for (int j = 0; j < M_; ++j) {
                    A_v(i, j, l) = An(i, j, l);
                    B_v(i, j, l) = Bn(i, j, l);
                    C_v(i, j, l) = Cn(i, j, l);
                }
            }
        }
    }

    // Final Stage: System is completely decoupled! Solve B'_l * X_l = D'_l
    for (int l = 0; l < nd; ++l) {
        double B_final[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
            CLOUDJ_UNROLL_4
            for (int j = 0; j < M_; ++j) {
                B_final[i][j] = B_v(i, j, l);
            }
        }
        double invB_final[M_][M_];
        invert_matrix_4x4_helper(B_final, invB_final);

        double D_final[M_];
        for (int i = 0; i < M_; ++i) D_final[i] = D_v(i, l);

        double X_final[M_];
        mat_vec_mult_4(invB_final, D_final, X_final);

        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
            rr(i, l) = X_final[i];
        }
    }

    // Extract boundary fluxes and populate output J-values
    fjtop = 0.0;
    fjbot = 0.0;
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
        fjtop += rr(i, 0) * WT[i];
        fjbot += rr(i, l_last) * WT[i];
    }

    for (int l = 0; l < nd; ++l) {
        double sum_rr = 0.0;
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
            sum_rr += rr(i, l) * WT[i];
        }
        fj(l, k_idx) = sum_rr;
    }
}
#endif

/**
 * @brief Main radiative transfer tridiagonal solver logic.
 * Translates subroutine BLKSLV in cldj_fjx_sub_mod.F90.
 */
inline void BLKSLV(
    mdspan_2d_mut fj,     // (N_, W_+W_r)
    mdspan_2d pomega,     // (M2_, N_) (from K-slice)
    mdspan_1d fz,         // (N_) (from K-slice)
    mdspan_1d ztau,       // (N_) (from K-slice)
    double fsbot,
    const std::array<double, 5>& rfl,
    const double pm[M_][M2_],
    const double pm0[M2_],
    double& fjtop,
    double& fjbot,
    std::array<double, 5>& fibot,
    int nd,
    int k_idx,            // Current wavelength index
    Workspace& ws         // Persistent pre-allocated workspace reference
) {
#if defined(CLOUDJ_USE_PCR)
    // Redirect cleanly to our Parallel Cyclic Reduction solver backend
    solve_pcr(fj, pomega, fz, ztau, fsbot, rfl, pm, pm0, fjtop, fjbot, fibot, nd, k_idx, ws);
    return;
#endif

    // Create mdspan wrappers directly mapping over persistent workspace buffers (zero allocation)
    mdspan_2d_mut a(ws.a_data.data(), M_, nd);
    mdspan_2d_mut c(ws.c_data.data(), M_, nd);
    mdspan_2d_mut h(ws.h_data.data(), M_, nd);
    mdspan_2d_mut rr(ws.rr_data.data(), M_, nd);

    mdspan_3d_mut b(ws.b_data.data(), M_, M_, nd);
    mdspan_3d_mut aa(ws.aa_data.data(), M_, M_, nd);
    mdspan_3d_mut cc(ws.cc_data.data(), M_, M_, nd);
    mdspan_3d_mut dd(ws.dd_data.data(), M_, M_, nd);

    // Generate block tri-diagonal system
    GEN_ID(pomega, fz, ztau, fsbot, rfl, pm, pm0, b, aa, cc, a, h, c, nd);

    // UPPER BOUNDARY L=1 (0-based: l=0)
    double E[M_][M_];
    for (int j = 0; j < M_; ++j) {
        for (int i = 0; i < M_; ++i) {
            E[i][j] = b(i, j, 0);
        }
    }

    solve_lu_4x4(E);

    for (int j = 0; j < M_; ++j) {
        for (int i = 0; i < M_; ++i) {
            dd(i, j, 0) = -E[i][0] * cc(0, j, 0) - E[i][1] * cc(1, j, 0) -
                          E[i][2] * cc(2, j, 0) - E[i][3] * cc(3, j, 0);
        }
        rr(j, 0) = E[j][0] * h(0, 0) + E[j][1] * h(1, 0) +
                   E[j][2] * h(2, 0) + E[j][3] * h(3, 0);
    }

    // CONTINUE THROUGH ALL DEPTH POINTS l=1 TO nd-2 (Fortran L=2 TO ND-1)
    for (int l = 1; l < nd - 1; ++l) {
        for (int j = 0; j < M_; ++j) {
            for (int i = 0; i < M_; ++i) {
                b(i, j, l) += a(i, l) * dd(i, j, l - 1);
            }
            h(j, l) -= a(j, l) * rr(j, l - 1);
        }

        for (int j = 0; j < M_; ++j) {
            for (int i = 0; i < M_; ++i) {
                E[i][j] = b(i, j, l);
            }
        }

        solve_lu_4x4(E);

        for (int j = 0; j < M_; ++j) {
            for (int i = 0; i < M_; ++i) {
                dd(i, j, l) = -E[i][j] * c(j, l);
            }
            rr(j, l) = E[j][0] * h(0, l) + E[j][1] * h(1, l) +
                       E[j][2] * h(2, l) + E[j][3] * h(3, l);
        }
    }

    // FINAL DEPTH POINT: l=nd-1 (Fortran L=ND)
    int l_last = nd - 1;
    for (int j = 0; j < M_; ++j) {
        for (int i = 0; i < M_; ++i) {
            b(i, j, l_last) += aa(i, 0, l_last) * dd(0, j, l_last - 1) +
                               aa(i, 1, l_last) * dd(1, j, l_last - 1) +
                               aa(i, 2, l_last) * dd(2, j, l_last - 1) +
                               aa(i, 3, l_last) * dd(3, j, l_last - 1);
        }
        h(j, l_last) -= aa(j, 0, l_last) * rr(0, l_last - 1) +
                        aa(j, 1, l_last) * rr(1, l_last - 1) +
                        aa(j, 2, l_last) * rr(2, l_last - 1) +
                        aa(j, 3, l_last) * rr(3, l_last - 1);
    }

    for (int j = 0; j < M_; ++j) {
        for (int i = 0; i < M_; ++i) {
            E[i][j] = b(i, j, l_last);
        }
    }

    solve_lu_4x4(E);

    for (int j = 0; j < M_; ++j) {
        rr(j, l_last) = E[j][0] * h(0, l_last) + E[j][1] * h(1, l_last) +
                        E[j][2] * h(2, l_last) + E[j][3] * h(3, l_last);
    }

    // BACK SOLUTION
    for (int l = nd - 2; l >= 0; --l) {
        for (int j = 0; j < M_; ++j) {
            rr(j, l) += dd(j, 0, l) * rr(0, l + 1) + dd(j, 1, l) * rr(1, l + 1) +
                        dd(j, 2, l) * rr(2, l + 1) + dd(j, 3, l) * rr(3, l + 1);
        }
    }

    // MEAN J & H (Fortran 1-based level structure)
    // 0-based L indices: L=0, 2, 4... are odd levels (Fortran 1, 3, 5...)
    // L=1, 3, 5... are even levels (Fortran 2, 4, 6...)
    for (int l = 0; l < nd; l += 2) {
        fj(l, k_idx) = rr(0, l) * WT[0] + rr(1, l) * WT[1] +
                       rr(2, l) * WT[2] + rr(3, l) * WT[3];
    }
    for (int l = 1; l < nd; l += 2) {
        fj(l, k_idx) = rr(0, l) * WT[0] * EMU[0] + rr(1, l) * WT[1] * EMU[1] +
                       rr(2, l) * WT[2] * EMU[2] + rr(3, l) * WT[3] * EMU[3];
    }

    // FJTOP diffuse flux out of top-of-atmosphere
    double sumt = rr(0, 0) * WT[0] * EMU[0] + rr(1, 0) * WT[1] * EMU[1] +
                  rr(2, 0) * WT[2] * EMU[2] + rr(3, 0) * WT[3] * EMU[3];
    fjtop = 4.0 * sumt;

    // Surface diffuse flux integration
    double sumb = rr(0, l_last) * WT[0] * EMU[0] + rr(1, l_last) * WT[1] * EMU[1] +
                  rr(2, l_last) * WT[2] * EMU[2] + rr(3, l_last) * WT[3] * EMU[3];
    
    double sumbr = rr(0, l_last) * WT[0] * EMU[0] * rfl[0] +
                   rr(1, l_last) * WT[1] * EMU[1] * rfl[1] +
                   rr(2, l_last) * WT[2] * EMU[2] * rfl[2] +
                   rr(3, l_last) * WT[3] * EMU[3] * rfl[3];

    double sumrf = WT[0] * EMU[0] * rfl[0] + WT[1] * EMU[1] * rfl[1] +
                   WT[2] * EMU[2] * rfl[2] + WT[3] * EMU[3] * rfl[3];

    double sumbx = (4.0 * sumbr + fsbot * rfl[4]) / (1.0 + 2.0 * sumrf);

    fjbot = 4.0 * sumb - sumbx;

    // FIBOT outputs: diffuse rays up/down
    fibot[4] = sumbx;
    for (int j = 0; j < 4; ++j) {
        fibot[j] = 2.0 * rr(j, l_last) - sumbx;
    }
}

} // namespace RadiativeSolver
} // namespace CloudJ

#endif // CLOUDJ_RADIATIVE_SOLVER_HPP
