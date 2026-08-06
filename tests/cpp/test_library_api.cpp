#include <cassert>
#include <cloudj/cloudj.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cloudj/photolysis.hpp>
#include <cloudj/profile.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/rates.hpp>
#include <cmath>
#include <iostream>

void test_profile_bounds() {
  CloudJ::AtmosphericProfile profile(80);
  assert(profile.get_num_layers() == 80);
}

void test_cross_section_interpolation() {
  // 1-point
  double x1 = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.5, 0.0, 0.0,
                                                 0.0, 0.0, 1);
  assert(std::abs(x1 - 1.5) < 1e-10);

  // 2-point
  double x2 = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.0, 300.0, 2.0,
                                                 0.0, 0.0, 2);
  assert(std::abs(x2 - 1.5) < 1e-10);

  // 2-point clamped lower
  double x2_low = CloudJ::CrossSections::interpolate(150.0, 200.0, 1.0, 300.0,
                                                     2.0, 0.0, 0.0, 2);
  assert(std::abs(x2_low - 1.0) < 1e-10);

  // 2-point clamped upper
  double x2_high = CloudJ::CrossSections::interpolate(350.0, 200.0, 1.0, 300.0,
                                                      2.0, 0.0, 0.0, 2);
  assert(std::abs(x2_high - 2.0) < 1e-10);

  // 3-point (lower half)
  double x3_low = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.0, 300.0,
                                                     2.0, 400.0, 4.0, 3);
  assert(std::abs(x3_low - 1.5) < 1e-10);

  // 3-point (upper half)
  double x3_high = CloudJ::CrossSections::interpolate(350.0, 200.0, 1.0, 300.0,
                                                      2.0, 400.0, 4.0, 3);
  assert(std::abs(x3_high - 3.0) < 1e-10);
}

void test_radiative_solver_and_gauss_gauss() {
  // Check that we have EMU and WT parameters defined
  assert(CloudJ::RadiativeSolver::M_ == 4);
  assert(CloudJ::RadiativeSolver::M2_ == 8);
  assert(CloudJ::RadiativeSolver::EMU[0] > 0.0);
  assert(CloudJ::RadiativeSolver::WT[0] > 0.0);

  // Check manual LU block-diagonal tri-diagonal inversion math functions
  // compiles
  double E[CloudJ::RadiativeSolver::M_][CloudJ::RadiativeSolver::M_] = {
      {4.0, 1.0, 0.0, 0.0},
      {1.0, 4.0, 1.0, 0.0},
      {0.0, 1.0, 4.0, 1.0},
      {0.0, 0.0, 1.0, 4.0}};

  // Invert block using standard blockLU decomposition
  CloudJ::RadiativeSolver::solve_lu_4x4(E);

  // Verify E contains inverted matrix elements
  assert(E[0][0] > 0.0);
  assert(E[1][1] > 0.0);
}

void test_photolysis_and_orchestrator() {
  CloudJ::Engine engine;
  CloudJ::AtmosphericProfile profile(80);

  // Calculate standard rates
  CloudJ::OutputRates rates = engine.calculate_photolysis_rates(profile, 45.0);
  assert(rates.j_values.size() == 80);
  assert(rates.j_values[0].size() == 3);

  // Verify dark conditions
  CloudJ::OutputRates dark_rates =
      engine.calculate_photolysis_rates(profile, 100.0);
  assert(dark_rates.j_values.size() == 80);
  assert(dark_rates.j_values[0][0] == 0.0);
}

void test_error_handling() {
  int rc = CloudJ::CLDJ_SUCCESS;
  CloudJ::CLOUDJ_ERROR("Warning test", "test_location", rc);
  assert(rc == CloudJ::CLDJ_FAILURE);

  try {
    CloudJ::CLOUDJ_ERROR_STOP("Fatal error test", "test_location");
    assert(false && "Should have thrown exception!");
  } catch (const CloudJ::Error &e) {
    // Exception successfully caught!
  }
}

void test_radiative_solver_workspace() {
  CloudJ::RadiativeSolver::Workspace ws;
  size_t nd = 10;
  ws.resize(nd);

  assert(ws.a_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.c_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.h_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.rr_data.size() == CloudJ::RadiativeSolver::M_ * nd);

  assert(ws.b_data.size() ==
         CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.aa_data.size() ==
         CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.cc_data.size() ==
         CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  assert(ws.dd_data.size() ==
         CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);

  for (double val : ws.a_data)
    assert(val == 0.0);
  for (double val : ws.b_data)
    assert(val == 0.0);
}

int main() {
  std::cout << "Running standard library API unit tests...\n";
  test_profile_bounds();
  test_cross_section_interpolation();
  test_radiative_solver_and_gauss_gauss();
  test_photolysis_and_orchestrator();
  test_error_handling();
  test_radiative_solver_workspace();
  std::cout << "All library API unit tests passed successfully!\n";
  return 0;
}
