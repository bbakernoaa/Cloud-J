#ifndef CLOUDJ_PHOTOLYSIS_HPP
#define CLOUDJ_PHOTOLYSIS_HPP

#include <cloudj/cross_sections.hpp>
#include <cmath>
#include <experimental/mdspan.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace CloudJ {
namespace Photolysis {

constexpr int W_ = 18; // standard wavelengths
constexpr int S_ = 27; // standard S-bins

struct SpecData {
  int nw = W_;
  int ns = S_;
  int njx = 0;
  std::vector<std::string> titlejx;
  std::vector<char> sqq; // 'p' or 't' interpolation variable
  std::vector<int> lqq;  // number of points (1, 2, or 3)
  std::vector<std::vector<double>>
      tqq; // interpolation temperatures/pressures per species [3]

  // cross-sections
  std::vector<std::vector<double>> qo2;              // [W_][3]
  std::vector<std::vector<double>> qo3;              // [W_][3]
  std::vector<std::vector<double>> q1d;              // [W_][3]
  std::vector<std::vector<std::vector<double>>> qqq; // [W_][3][species]

  // Pre-computed reciprocal temperature span intervals
  std::vector<double>
      inv_t12; // [species] reciprocal of (tqq[species][1] - tqq[species][0])
  std::vector<double>
      inv_t23; // [species] reciprocal of (tqq[species][2] - tqq[species][1])
};

/**
 * @brief Interpolates and accumulates photolysis rates (J-values) for a column.
 * Matches JRATET in cldj_fjx_sub_mod.F90.
 */
inline void
JRATET(const std::vector<double> &ppj, // pressure edges [lu + 1]
       const std::vector<double> &ttj, // mid-layer temperatures [lu + 1]
       std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>,
                                 std::experimental::layout_left>
           fff,                                 // mean actinic fluxes [W_][lu]
       std::vector<std::vector<double>> &valjl, // [lu][njx]
       const SpecData &spec, int lu, int njxu) {
  valjl.assign(lu, std::vector<double>(njxu, 0.0));

  for (int l = 0; l < lu; ++l) {
    double tt = ttj[l];
    double pp;
    if (l == 0) {
      pp = ppj[0];
    } else {
      pp = (ppj[l] + ppj[l + 1]) * 0.5;
    }

    // zero bin-11 below 100 hPa matching the O2 e-fold limit
    if (pp > 100.0) {
      fff(10, l) = 0.0; // 0-based index 10 corresponds to bin 11
    }

    std::vector<double> valj(spec.njx, 0.0);

    // Calculate O2, O3, and O3(1D) photolysis rates (reactions 0, 1, 2)
    for (int k = 0; k < W_; ++k) {
      double qo2tot = CrossSections::interpolate(
          tt, spec.tqq[0][0], spec.qo2[k][0], spec.tqq[0][1], spec.qo2[k][1],
          spec.tqq[0][2], spec.qo2[k][2], spec.lqq[0], spec.inv_t12[0],
          spec.inv_t23[0]);

      double qo3tot = CrossSections::interpolate(
          tt, spec.tqq[1][0], spec.qo3[k][0], spec.tqq[1][1], spec.qo3[k][1],
          spec.tqq[1][2], spec.qo3[k][2], spec.lqq[1], spec.inv_t12[1],
          spec.inv_t23[1]);

      double qo31dy = CrossSections::interpolate(
          tt, spec.tqq[2][0], spec.q1d[k][0], spec.tqq[2][1], spec.q1d[k][1],
          spec.tqq[2][2], spec.q1d[k][2], spec.lqq[2], spec.inv_t12[2],
          spec.inv_t23[2]);

      double qo31d = qo31dy * qo3tot;

      valj[0] += qo2tot * fff(k, l);
      valj[1] += qo3tot * fff(k, l);
      valj[2] += qo31d * fff(k, l);
    }

    // Calculate photolysis rates for reactions 4 to NJX (indices 3 to NJX-1)
    for (int j = 3; j < spec.njx; ++j) {
      double var = (spec.sqq[j] == 'p') ? pp : tt;
      for (int k = 0; k < W_; ++k) {
        double qqqt = CrossSections::interpolate(
            var, spec.tqq[j][0], spec.qqq[k][0][j], spec.tqq[j][1],
            spec.qqq[k][1][j], spec.tqq[j][2], spec.qqq[k][2][j], spec.lqq[j],
            spec.inv_t12[j], spec.inv_t23[j]);
        valj[j] += qqqt * fff(k, l);
      }
    }

    for (int j = 0; j < spec.njx; ++j) {
      valjl[l][j] = valj[j];
    }
  }
}

} // namespace Photolysis
} // namespace CloudJ

#endif // CLOUDJ_PHOTOLYSIS_HPP
