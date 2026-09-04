#ifndef CLOUDJ_CLOUDJ_HPP
#define CLOUDJ_CLOUDJ_HPP

#include <array>
#include <cloudj/cloud_jx.hpp>
#include <cloudj/context.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cloudj/init.hpp>
#include <cloudj/kokkos_backend.hpp>
#include <cloudj/osa.hpp>
#include <cloudj/photo_jx.hpp>
#include <cloudj/photolysis.hpp>
#include <cloudj/profile.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/rates.hpp>
#include <cloudj/state.hpp>
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
  CloudJState state_;
  Photolysis::SpecData spec_data;
  RadiativeSolver::Workspace solver_ws;
  bool initialized_ = false;

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
    // (flat layout: tqq[j * 3 + t], qo2/qo3/q1d[k * 3 + t])
    spec_data.tqq.resize(3 * 3);
    for (int j = 0; j < 3; ++j) {
      spec_data.tqq[j * 3 + 0] = 200.0;
      spec_data.tqq[j * 3 + 1] = 300.0;
      spec_data.tqq[j * 3 + 2] = 400.0;
    }

    spec_data.qo2.resize(Photolysis::W_ * 3);
    spec_data.qo3.resize(Photolysis::W_ * 3);
    spec_data.q1d.resize(Photolysis::W_ * 3);
    for (int k = 0; k < Photolysis::W_; ++k) {
      spec_data.qo2[k * 3 + 0] = 1e-20;
      spec_data.qo2[k * 3 + 1] = 2e-20;
      spec_data.qo2[k * 3 + 2] = 3e-20;

      spec_data.qo3[k * 3 + 0] = 1e-19;
      spec_data.qo3[k * 3 + 1] = 2e-19;
      spec_data.qo3[k * 3 + 2] = 3e-19;

      spec_data.q1d[k * 3 + 0] = 0.1;
      spec_data.q1d[k * 3 + 1] = 0.5;
      spec_data.q1d[k * 3 + 2] = 0.9;
    }

    // Setup pre-computed reciprocal temperature span intervals
    spec_data.inv_t12.assign(3, 0.0);
    spec_data.inv_t23.assign(3, 0.0);
    for (size_t j = 0; j < 3; ++j) {
      spec_data.inv_t12[j] =
          1.0 / (spec_data.tqq[j * 3 + 1] - spec_data.tqq[j * 3 + 0]);
      spec_data.inv_t23[j] =
          1.0 / (spec_data.tqq[j * 3 + 2] - spec_data.tqq[j * 3 + 1]);
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
    for (int l = 0; l < l1u; ++l) {
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
    // sections. Passing the solved mean actinic flux (fjact) to evaluate final
    // J-values. JRATET now writes a flat row-major buffer [l*njx + j]; unpack
    // it into the public OutputRates::j_values nested layout.
    std::vector<double> valjl_flat;
    Photolysis::JRATET(ppj, ttj, fjact, valjl_flat, spec_data, lu,
                       spec_data.njx);
    for (int l = 0; l < lu; ++l) {
      for (int j = 0; j < spec_data.njx; ++j) {
        rates.j_values[l][j] = valjl_flat[l * spec_data.njx + j];
      }
    }

    return rates;
  }

  /**
   * @brief Initialize the Cloud-J engine by loading all lookup tables into
   * CloudJState. Wraps Init::INIT_CLDJ().
   *
   * @param amiroot    True if this is the root process (controls diagnostic output)
   * @param nlevels    Number of CTM levels
   * @param nlevels_with_cloud Number of CTM levels that include cloud data
   * @param njxu       Number of photolysis reactions requested by the caller
   * @param atau       Geometric factor for inserted sub-layers (default 1.120)
   * @param atau0      Minimum OD threshold for inserted layers (default 0.005)
   * @param nwbin      Number of wavelength bins to use (8, 12, or 18)
   * @param cldflag    Cloud overlap scheme flag (1-8, not 4)
   * @param cldcor     Cloud decorrelation parameter
   * @param lnrg       Max-overlap group selection (0, 3, or 6)
   * @param atm0       Atmosphere type selector for reference profiles
   * @param use_h2o_uv_abs Enable H2O UV absorption
   * @param njxx       Output: actual number of photolysis reactions available
   * @param datadir    Optional path to a directory containing on-disk override
   *                   tables (FJX_spec.dat, FJX_scat-cld.dat, etc. -- see
   *                   init.hpp for the full filename list). If empty (default)
   *                   or a specific file is not found there, falls back to the
   *                   embedded compiled-in tables. A message is printed to
   *                   stderr only when a disk file is actually used.
   * @return           CLDJ_SUCCESS (0) on success, non-zero on error
   */
  int init(bool amiroot, int nlevels, int nlevels_with_cloud,
           int njxu, double atau, double atau0, int nwbin,
           int cldflag, double cldcor, int lnrg, int atm0,
           bool use_h2o_uv_abs, int& njxx,
           const std::string& datadir = "") {
    int rc = CLDJ_SUCCESS;
    std::vector<std::string> titlejxx;

    Init::INIT_CLDJ(
        amiroot, datadir,
        nlevels, nlevels_with_cloud,
        titlejxx, njxu,
        atau, atau0, nwbin,
        cldflag, cldcor, lnrg, atm0,
        use_h2o_uv_abs,
        njxx, state_, rc);

    if (rc == CLDJ_SUCCESS) {
      initialized_ = true;
    }
    return rc;
  }

  /**
   * @brief Execute the full CLOUD_JX pipeline for a single atmospheric column.
   * Wraps the global CLOUD_JX() function using the Engine's internal CloudJState.
   *
   * @param u0         cos(solar zenith angle)
   * @param sza        Solar zenith angle (degrees)
   * @param rfl_flat   Surface reflectivity array [5 * (W_ + W_r)]
   * @param solf       Solar flux factor
   * @param lprtj      Enable diagnostic printing
   * @param ppp        Pressure edges [L1U+1]
   * @param zzz        Height edges [L1U+1] (cm)
   * @param ttt        Temperature per layer [L1U]
   * @param hhh        H2O column per layer [L1U]
   * @param ddd        Air density per layer [L1U]
   * @param rrr        Relative humidity per layer [L1U]
   * @param ooo        O3 column per layer [L1U]
   * @param ccc        CH4 column per layer [L1U]
   * @param lwp        Liquid water path per layer [L1U]
   * @param iwp        Ice water path per layer [L1U]
   * @param reffl      Liquid cloud effective radius [L1U]
   * @param reffi      Ice cloud effective radius [L1U]
   * @param cldf       Cloud fraction per layer [L1U]
   * @param cldiw      Cloud type index per layer [L1U]
   * @param cldcor_in  Cloud decorrelation parameter
   * @param aersp      Aerosol species paths [AN_ * L1U]
   * @param ndxaer     Aerosol type indices [AN_ * L1U]
   * @param l1u        Number of levels (L1U = L_ + 1)
   * @param anu        Number of aerosol types in layer (AN_)
   * @param njxu       Number of photolysis reactions
   * @param valjxx     Output: J-values [(L1U-1) * NJXU]
   * @param skperd     Output: heating rates [(S_+2) * L1U]
   * @param swmsq      Output: 6 solar flux diagnostics
   * @param od18       Output: per-layer OD at 18 bins [L1U]
   * @param iran       Random number seed index
   * @param nica       Output: number of ICAs computed
   * @param jcount     Output: number of PHOTO_JX calls made
   * @param ldark      Output: true if column is dark (SZA > limit)
   * @param wtqca      Output: quadrature weights [NQD_]
   * @param rc         Output: return code (0=success)
   * @param dir_sfc_flux   Optional output: direct surface flux [W_]
   * @param diff_sfc_flux  Optional output: diffuse surface flux [W_]
   * @param dep_flux       Optional output: deposition flux [W_]
   * @param diff_top_flux  Optional output: diffuse top-of-atmosphere flux [W_]
   */
  void cloud_jx(
      double u0, double sza,
      const double* rfl_flat, double solf,
      bool lprtj,
      const double* ppp, const double* zzz,
      const double* ttt, const double* hhh,
      const double* ddd, const double* rrr,
      const double* ooo, const double* ccc,
      const double* lwp, const double* iwp,
      const double* reffl, const double* reffi,
      const double* cldf, const int* cldiw,
      double cldcor_in,
      const double* aersp, const int* ndxaer,
      int l1u, int anu, int njxu,
      double* valjxx, double* skperd, double* swmsq, double* od18,
      int iran, int& nica, int& jcount, bool& ldark,
      double* wtqca,
      int& rc,
      double* dir_sfc_flux = nullptr,
      double* diff_sfc_flux = nullptr,
      double* dep_flux = nullptr,
      double* diff_top_flux = nullptr) {

    CLOUD_JX(u0, sza, rfl_flat, solf, lprtj,
             ppp, zzz, ttt, hhh, ddd, rrr, ooo, ccc,
             lwp, iwp, reffl, reffi, cldf, cldiw,
             cldcor_in,
             aersp, ndxaer,
             l1u, anu, njxu,
             valjxx, skperd, swmsq, od18,
             iran, nica, jcount, ldark,
             wtqca,
             state_, rc,
             dir_sfc_flux, diff_sfc_flux, dep_flux, diff_top_flux);
  }

  /**
   * @brief Returns a const reference to the internal CloudJState.
   */
  const CloudJState& get_state() const noexcept { return state_; }
};

} // namespace CloudJ

#endif // CLOUDJ_CLOUDJ_HPP
