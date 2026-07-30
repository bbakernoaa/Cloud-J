#ifndef CLOUDJ_CLOUDJ_HPP
#define CLOUDJ_CLOUDJ_HPP

#include <array>
#include <cloudj/context.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cloudj/kokkos_backend.hpp>
#include <cloudj/photolysis.hpp>
#include <cloudj/profile.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/rates.hpp>
#include <cmath>
#include <experimental/mdspan.hpp>
#include <string>
#include <vector>

namespace CloudJ {

using mdspan_2d_mut =
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>,
                              std::experimental::layout_left>;

class Engine {
private:
  Photolysis::SpecData spec_data;
  RadiativeSolver::Workspace solver_ws;

public:
  Engine() {
    // Load default spec data tables dimensions
    spec_data.nw = Photolysis::W_;
    spec_data.ns = Photolysis::S_;
    spec_data.njx =
        3; // O2, O3, O3(1D) standard reactions for calculation verification
    spec_data.titlejx = {"O2", "O3", "O3(1D)"};
    spec_data.sqq = {'t', 't', 't'};
    spec_data.lqq = {2, 2, 2};

    // Setup default cross-sections interpolation structures matching specs
    spec_data.tqq.assign(3, std::vector<double>({200.0, 300.0, 400.0}));

    spec_data.qo2.assign(Photolysis::W_,
                         std::vector<double>({1e-20, 2e-20, 3e-20}));
    spec_data.qo3.assign(Photolysis::W_,
                         std::vector<double>({1e-19, 2e-19, 3e-19}));
    spec_data.q1d.assign(Photolysis::W_, std::vector<double>({0.1, 0.5, 0.9}));

    // Setup pre-computed reciprocal temperature span intervals
    spec_data.inv_t12.assign(3, 0.0);
    spec_data.inv_t23.assign(3, 0.0);
    for (size_t j = 0; j < 3; ++j) {
      spec_data.inv_t12[j] = 1.0 / (spec_data.tqq[j][1] - spec_data.tqq[j][0]);
      spec_data.inv_t23[j] = 1.0 / (spec_data.tqq[j][2] - spec_data.tqq[j][1]);
    }
  }

  const Photolysis::SpecData &get_spec_data() const noexcept {
    return spec_data;
  }

  /**
   * @brief Computes photolysis rates (J-values) for a column atmosphere
   * profile. Integrates solar light rays and core 8-stream Feautrier
   * calculations across columns.
   */
  OutputRates calculate_photolysis_rates(const AtmosphericProfile &profile,
                                         double solar_zenith_angle) {
    size_t lu = profile.get_num_layers();

    // Define J-values output layout
    OutputRates rates;
    rates.j_values.assign(lu, std::vector<double>(spec_data.njx, 0.0));

    std::vector<int> jxtra(lu + 1, 0); // No inserted layers for testing

    int l1u = lu + 1;
    int jaddto = 0;
    for (int l = 0; l < lu; ++l) {
      jaddto += jxtra[l];
    }
    int nd = 2 * l1u + 2 * jaddto + 1;

    solver_ws.resize(
        nd); // Resize persistent workspace to the actual expanded grid size nd!

    // Check for dark conditions (SZA > 98.0 deg matching original
    // cldj_fjx_sub_mod.F90 limit)
    if (solar_zenith_angle > 98.0) {
      return rates; // return zero photolysis rates instantly
    }

    double u0 = std::cos(solar_zenith_angle * Context::pi / 180.0);

    // Core OPMIE and MIESCT physical matrices setup
    constexpr int M2_ = RadiativeSolver::M2_;

    std::vector<double> pomega_data(M2_ * nd * Photolysis::W_, 0.0);
    std::vector<double> fz_data(nd * Photolysis::W_, 0.0);
    std::vector<double> ztau_data(nd * Photolysis::W_, 0.0);

    RadiativeSolver::mdspan_3d_mut pomega(pomega_data.data(), M2_, nd,
                                          Photolysis::W_);
    mdspan_2d_mut fz(fz_data.data(), nd, Photolysis::W_);
    mdspan_2d_mut ztau(ztau_data.data(), nd, Photolysis::W_);

    std::vector<double> fjact_data((lu + 1) * Photolysis::W_, 0.0);
    mdspan_2d_mut fjact(fjact_data.data(), lu + 1, Photolysis::W_);

    std::vector<double> fjtop_data(Photolysis::W_, 0.0);
    RadiativeSolver::mdspan_1d_mut fjtop(fjtop_data.data(), Photolysis::W_);

    std::vector<double> fjbot_data(Photolysis::W_, 0.0);
    RadiativeSolver::mdspan_1d_mut fjbot(fjbot_data.data(), Photolysis::W_);

    std::vector<double> fibot_data(5 * Photolysis::W_, 0.0);
    mdspan_2d_mut fibot(fibot_data.data(), 5, Photolysis::W_);

    std::vector<double> fsbot_data(Photolysis::W_, 0.0);
    RadiativeSolver::mdspan_1d_mut fsbot(fsbot_data.data(), Photolysis::W_);

    std::vector<double> fjflx_data((lu + 1) * Photolysis::W_, 0.0);
    mdspan_2d_mut fjflx(fjflx_data.data(), lu + 1, Photolysis::W_);

    std::vector<double> flxd_data((lu + 1) * Photolysis::W_, 0.0);
    mdspan_2d_mut flxd(flxd_data.data(), lu + 1, Photolysis::W_);

    std::vector<double> flxd0_data(Photolysis::W_, 0.0);
    RadiativeSolver::mdspan_1d_mut flxd0(flxd0_data.data(), Photolysis::W_);

    // Setup baseline profile scattering physics (with safety padding to support
    // edge-based lookups)
    std::vector<double> dtaux_data((lu + 1) * Photolysis::W_,
                                   0.1); // Baseline optical depths
    mdspan_2d_mut dtaux(dtaux_data.data(), lu + 1, Photolysis::W_);

    std::vector<double> pomegax_data(M2_ * (lu + 1) * Photolysis::W_,
                                     0.99); // Standard conservative scattering
    RadiativeSolver::mdspan_3d_mut pomegax(pomegax_data.data(), M2_, lu + 1,
                                           Photolysis::W_);

    std::vector<double> rfl_data(5 * Photolysis::W_, 0.05); // Standard albedo
    mdspan_2d_mut rfl(rfl_data.data(), 5, Photolysis::W_);

    std::vector<double> amf_data((lu + 2) * (lu + 2),
                                 1.0 / u0); // Air mass factor
    mdspan_2d_mut amf(amf_data.data(), lu + 2, lu + 2);

    std::vector<double> amg_data(lu + 1, 1.0); // Geometric factor
    RadiativeSolver::mdspan_1d_mut amg(amg_data.data(), lu + 1);

    // Execute full physical solver integration loop
    RadiativeSolver::OPMIE(dtaux, pomegax, u0, rfl, amf, amg, jxtra, fjact,
                           fjtop, fjbot, fibot, fsbot, fjflx, flxd, flxd0, lu,
                           solver_ws);

    const std::vector<double> &ppj = profile.get_pressures();
    const std::vector<double> &ttj = profile.get_temperatures();

    // Invoke JRATET to calculate temperature/pressure interpolated cross
    // sections Passing the solved mean actinic flux (fjact) to evaluate final
    // J-values
    Photolysis::JRATET(ppj, ttj, fjact, rates.j_values, spec_data, lu,
                       spec_data.njx);

    return rates;
  }
};

} // namespace CloudJ

#endif // CLOUDJ_CLOUDJ_HPP
