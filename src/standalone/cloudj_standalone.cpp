// Cloud-J v8.0 C++ Standalone Driver
// Mirrors the Fortran cldj_standalone.F90 to produce identical output.

#ifndef MODEL_STANDALONE
#define MODEL_STANDALONE
#endif

#include <algorithm>
#include <chrono>
#include <cloudj/cloudj.hpp>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Convert a double to Fortran's e9.2 scientific formatting cleanly
std::string format_fortran_e9_2(double val) {
  if (std::abs(val) < 1e-99) {
    return " 0.00E+00";
  }
  std::ostringstream ss;
  ss << std::scientific << std::uppercase << std::setprecision(2) << val;
  std::string s = ss.str();

  size_t e_pos = s.find('E');
  if (e_pos == std::string::npos) {
    return " 0.00E+00";
  }

  std::string mantissa = s.substr(0, e_pos);
  std::string exp_part = s.substr(e_pos + 1);

  char exp_sign = '+';
  if (exp_part[0] == '-' || exp_part[0] == '+') {
    exp_sign = exp_part[0];
    exp_part = exp_part.substr(1);
  }

  int exp_val = std::stoi(exp_part);
  std::ostringstream exp_ss;
  exp_ss << exp_sign << std::setw(2) << std::setfill('0') << exp_val;

  std::string result = mantissa + "E" + exp_ss.str();
  if (result[0] != '-') {
    result = " " + result;
  }
  while (result.length() < 9) {
    result = " " + result;
  }
  return result;
}

int main(int argc, char *argv[]) {
  std::string tables_dir = "tables";
  std::string output_path = "cpp_actual_output.txt";
  // bench_iters: number of in-memory repeated calls to engine.cloud_jx() to
  // time, excluding one-time init/atmosphere-setup cost. 0 = disabled
  // (default single-pass output-comparison mode).
  int bench_iters = 0;
  // init_tables_dir: optional directory of on-disk override tables passed to
  // engine.init() for the engine's internal spectral/aerosol/cloud tables
  // (FJX_spec.dat, FJX_scat-cld.dat, etc). Empty by default = embedded
  // (compiled-in) tables, matching current behavior. This is distinct from
  // `tables_dir` above, which locates the atmos_PTClds.dat profile file that
  // this driver reads directly and is unrelated to the engine's tables.
  std::string init_tables_dir = "";

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--output" && i + 1 < argc) {
      output_path = argv[++i];
    } else if (arg == "--bench-iters" && i + 1 < argc) {
      bench_iters = std::atoi(argv[++i]);
    } else if (arg == "--benchmark") {
      // Backward-compatible alias for benchmark/run_benchmark.py's
      // --benchmark-iters harness flag, which forwards a bare "--benchmark"
      // to this executable. Use a small default iteration count so the
      // harness still measures the corrected cloud_jx() throughput instead
      // of being a no-op.
      if (bench_iters == 0) bench_iters = 1000;
    } else if (arg == "--tables-dir" && i + 1 < argc) {
      // Only configures the engine's internal table loading (init_tables_dir).
      // We deliberately do NOT also repoint the existing `tables_dir` variable
      // (used for reading atmos_PTClds.dat) here: that file is not part of the
      // 10-table mapping this flag documents, and changing its resolution
      // based on this flag would be a surprising, undocumented side effect
      // for a flag named "tables-dir" in the context of engine table loading.
      init_tables_dir = argv[++i];
    }
  }

  // =====================================================================
  // Constants from state.hpp (MODEL_STANDALONE)
  // =====================================================================
  constexpr int L_  = CloudJ::CloudJState::L_;   // 57
  constexpr int L1_ = CloudJ::CloudJState::L1_;  // 58
  constexpr int L2_ = CloudJ::CloudJState::L2_;  // 59
  constexpr int LWEPAR = CloudJ::CloudJState::LWEPAR; // 34

  // =====================================================================
  // Engine initialization (mirrors Fortran INIT_CLDJ)
  // =====================================================================
  CloudJ::Engine engine;

  double ATAU_in    = 1.050;
  double ATAU0_in   = 0.005;
  double CLDCOR_in  = 0.33;
  int    NWBIN_in   = 18;
  int    LNRG_in    = 6;
  int    ATM0_in    = 1;
  int    CLDFLAG_in = 7;
  bool   use_h2o_uv_abs = true;

  int NJXX = 0;
  int rc = engine.init(
      true,       // amiroot
      L_,         // nlevels = 57
      LWEPAR,     // nlevels_with_cloud = 34
      CloudJ::JVN_,  // njxu
      ATAU_in, ATAU0_in, NWBIN_in,
      CLDFLAG_in, CLDCOR_in, LNRG_in, ATM0_in,
      use_h2o_uv_abs, NJXX,
      init_tables_dir);  // empty string by default = embedded tables

  if (rc != 0) {
    std::cerr << "Error: INIT_CLDJ failed with rc=" << rc << "\n";
    return 1;
  }

  // =====================================================================
  // Read atmosphere file: tables/atmos_PTClds.dat
  // =====================================================================
  std::string profile_filename = tables_dir + "/atmos_PTClds.dat";
  std::ifstream infile(profile_filename);
  if (!infile.is_open()) {
    profile_filename = "./tables/atmos_PTClds.dat";
    infile.open(profile_filename);
    if (!infile.is_open()) {
      std::cerr << "Error: Could not open profile file tables/atmos_PTClds.dat\n";
      return 1;
    }
  }

  std::string line;
  int MONTH = 0, ILAT = 0;
  double PSURF = 0.0;
  double ALBEDO[5] = {};
  double WIND = 0.0, CHLR = 0.0;

  // Line 1: title (skip)
  std::getline(infile, line);
  // Line 2: MONTH, ILAT (format 2i5)
  std::getline(infile, line);
  MONTH = std::stoi(line.substr(0, 5));
  ILAT  = std::stoi(line.substr(5, 5));
  // Line 3: PSURF (format f5.0)
  std::getline(infile, line);
  PSURF = std::stod(line.substr(0, 5));
  // Line 4: ALBEDO(5) (format f5.2)
  std::getline(infile, line);
  ALBEDO[4] = std::stod(line.substr(0, 5));
  // Line 5: ALBEDO(1:4) (format 4f5.2)
  std::getline(infile, line);
  for (int i = 0; i < 4; ++i)
    ALBEDO[i] = std::stod(line.substr(i * 5, 5));
  // Line 6: WIND, CHLR
  std::getline(infile, line);
  {
    std::istringstream ss(line);
    ss >> WIND >> CHLR;
  }
  // Line 7: header (skip)
  std::getline(infile, line);

  // Read L1_ (58) levels of atmosphere data
  double ETAA[L2_] = {}, ETAB[L2_] = {};
  double TINP[L1_] = {}, RHINP[L1_] = {}, ZOFL[L1_] = {};
  double AER1[L1_] = {}, AER2[L1_] = {};
  int    NAA1[L1_] = {}, NAA2[L1_] = {};

  for (int L = 0; L < L1_; ++L) {
    std::getline(infile, line);
    std::istringstream ss(line);
    int idx;
    ss >> idx >> ETAA[L] >> ETAB[L] >> TINP[L] >> RHINP[L]
       >> ZOFL[L] >> AER1[L] >> NAA1[L] >> AER2[L] >> NAA2[L];
  }

  // Read cloud header
  std::getline(infile, line);

  // Read LWEPAR cloud layers (reversed: LWEPAR down to 1)
  double CLDFRW[LWEPAR] = {}, CLDLWCW[LWEPAR] = {}, CLDIWCW[LWEPAR] = {};
  for (int L = LWEPAR - 1; L >= 0; --L) {
    std::getline(infile, line);
    std::istringstream ss(line);
    int idx;
    ss >> idx >> CLDFRW[L] >> CLDLWCW[L] >> CLDIWCW[L];
  }
  infile.close();

  // =====================================================================
  // Compute pressure edges
  // =====================================================================
  ETAA[L2_ - 1] = 0.0;
  ETAB[L2_ - 1] = 0.0;
  double PPP[L2_] = {};
  for (int L = 0; L < L2_; ++L) {
    PPP[L] = ETAA[L] + ETAB[L] * PSURF;
  }

  // =====================================================================
  // Call ACLIM_FJX for O3, T, CH4 climatologies
  // =====================================================================
  double TTT[L1_] = {}, OOO[L1_] = {}, CH4[L1_] = {};
  double O3[L1_] = {};

  double YLAT = static_cast<double>(ILAT);
  CloudJ::PhotoJX::ACLIM_FJX(MONTH, YLAT, PPP, TTT, O3, CH4, L1_, engine.get_state());

  // Override T and RH from file (keep O3 from climatology)
  double RRR[L1_] = {};
  for (int L = 0; L < L1_; ++L) {
    TTT[L] = TINP[L];
    RRR[L] = RHINP[L];
  }

  // =====================================================================
  // Compute altitudes, densities, O3/CH4 columns
  // =====================================================================
  double ZZZ[L2_] = {};
  double DDD[L1_] = {};
  double CCC[L1_] = {};

  ZZZ[0] = 16.0e5 * std::log10(1013.25 / PPP[0]);  // cm

  for (int L = 0; L < L_; ++L) {
    DDD[L] = (PPP[L] - PPP[L + 1]) * CloudJ::MASFAC;
    double SCALEH = 1.3806e-19 * CloudJ::MASFAC * TTT[L];
    ZZZ[L + 1] = ZZZ[L] - (std::log(PPP[L + 1] / PPP[L]) * SCALEH);
    OOO[L] = DDD[L] * O3[L] * 1.0e-6;
    CCC[L] = DDD[L] * CH4[L] * 1.0e-9;
  }
  // Top layer L_ (index L_ = L1_-1 in 0-based)
  {
    int L = L_;  // = 57, the extra top layer
    ZZZ[L + 1] = ZZZ[L] + CloudJ::ZZHT;
    DDD[L] = (PPP[L] - PPP[L + 1]) * CloudJ::MASFAC;
    OOO[L] = DDD[L] * O3[L] * 1.0e-6;
    CCC[L] = DDD[L] * CH4[L] * 1.0e-9;
  }

  // =====================================================================
  // H2O profile
  // =====================================================================
  double HHH[L1_] = {};
  double HHH0 = 0.030;
  for (int L = 0; L < L1_; ++L) {
    HHH[L] = DDD[L] * std::max(HHH0 * std::exp(-ZZZ[L] / 2.2e5), 2.0e-6);
  }

  // =====================================================================
  // Aerosols setup
  // Layout: column-major [L1_][AN_], accessed as arr[L + L1_ * M]
  // =====================================================================
  double AERSP[L1_ * CloudJ::AN_] = {};
  int    NDXAER[L1_ * CloudJ::AN_] = {};

  for (int L = 0; L < L_; ++L) {
    NDXAER[L + L1_ * 0] = NAA1[L];
    AERSP[L + L1_ * 0]  = AER1[L];
    NDXAER[L + L1_ * 1] = NAA2[L];
    AERSP[L + L1_ * 1]  = AER2[L];
  }

  // =====================================================================
  // Cloud processing
  // =====================================================================
  double CLF[L1_] = {};
  double WLC[L_] = {}, WIC[L_] = {};
  double LWP[L1_] = {}, IWP[L1_] = {};
  double REFFL[L1_] = {}, REFFI[L1_] = {};
  int    CLDIW[L1_] = {};
  int    LTOP = LWEPAR;

  // Check if all cloud fractions are negligible
  double max_clf = 0.0;
  for (int L = 0; L < LWEPAR; ++L)
    max_clf = std::max(max_clf, CLDFRW[L]);

  if (max_clf <= 0.005) {
    // No clouds
    std::memset(IWP, 0, sizeof(IWP));
    std::memset(REFFI, 0, sizeof(REFFI));
    std::memset(LWP, 0, sizeof(LWP));
    std::memset(REFFL, 0, sizeof(REFFL));
  }

  for (int L = 0; L < LTOP; ++L) {
    CLDIW[L] = 0;
    double CF = CLDFRW[L];
    if (CF > 0.005) {
      CLF[L] = CF;
      WLC[L] = CLDLWCW[L] / CF;
      WIC[L] = CLDIWCW[L] / CF;
      if (WLC[L] > 1.0e-11) CLDIW[L] = 1;
      if (WIC[L] > 1.0e-11) CLDIW[L] = CLDIW[L] + 2;
    } else {
      CLF[L] = 0.0;
      WLC[L] = 0.0;
      WIC[L] = 0.0;
    }
  }

  // Derive R-effective for clouds
  for (int L = 0; L < LTOP; ++L) {
    // Ice clouds
    if (WIC[L] > 1.0e-12) {
      double PDEL = PPP[L] - PPP[L + 1];
      double ZDEL = (ZZZ[L + 1] - ZZZ[L]) * 0.01;  // m
      IWP[L] = 1000.0 * WIC[L] * PDEL * CloudJ::G100;  // g/m2
      double ICWC = IWP[L] / ZDEL;  // g/m3
      REFFI[L] = 164.0 * std::pow(ICWC, 0.23);
    } else {
      IWP[L] = 0.0;
      REFFI[L] = 0.0;
    }
    // Water clouds
    if (WLC[L] > 1.0e-12) {
      double PMID = 0.5 * (PPP[L] + PPP[L + 1]);
      double PDEL = PPP[L] - PPP[L + 1];
      double F1 = 0.005 * (PMID - 610.0);
      F1 = std::min(1.0, std::max(0.0, F1));
      LWP[L] = 1000.0 * WLC[L] * PDEL * CloudJ::G100;  // g/m2
      REFFL[L] = 9.6 * F1 + 12.68 * (1.0 - F1);
    } else {
      LWP[L] = 0.0;
      REFFL[L] = 0.0;
    }
  }

  // Copy CLF back to CLDFRW for use in SZA loop
  for (int L = 0; L < LTOP; ++L) {
    CLDFRW[L] = CLF[L];
  }

  // Total spectral bins
  constexpr int WW = CloudJ::W_ + CloudJ::W_r;  // = 18 + 0 = 18

  // =====================================================================
  // In-memory benchmark mode: exercises the CORRECT, numerically-verified
  // engine.cloud_jx() path repeatedly, excluding the one-time init and
  // atmosphere-setup cost above (which a host model would only pay once
  // per column/timestep before calling cloud_jx). Reuses the same output
  // buffers across iterations, matching how a host model reuses buffers
  // across timesteps.
  // =====================================================================
  if (bench_iters > 0) {
    const CloudJ::CloudJState& state = engine.get_state();

    int NSZA = 30;  // mid-range, representative SZA
    double SZA = static_cast<double>(NSZA);
    int IRAN = 1;
    double SOLF = 1.0;
    double U0 = std::cos(SZA * CloudJ::CPI180);

    double ANGLES[5];
    ANGLES[0] = CloudJ::EMU[0];
    ANGLES[1] = CloudJ::EMU[1];
    ANGLES[2] = CloudJ::EMU[2];
    ANGLES[3] = CloudJ::EMU[3];
    ANGLES[4] = U0;

    double RFL[5 * WW] = {};
    for (int K = 0; K < CloudJ::NS2; ++K) {
      double WAVEL = state.WL[K];
      double OSA_dir[5] = {};
      CloudJ::OSA::FJX_OSA(WAVEL, WIND, CHLR, ANGLES, OSA_dir);
      for (int J = 0; J < 5; ++J) {
        RFL[J + 5 * K] = OSA_dir[J];
        // Override OSA with read-in ALBEDO values (same as normal mode)
        RFL[J + 5 * K] = ALBEDO[J];
      }
    }

    // Disable diagnostic printing so I/O does not contaminate timing.
    bool LPRTJ = false;

    // Output buffers, reused (overwritten) across all iterations.
    double VALJXX[L_ * CloudJ::JVN_] = {};
    double SKPERD[(CloudJ::S_ + 2) * L1_] = {};
    double SWMSQ[6] = {};
    double OD18[L1_] = {};
    double WTQCA[CloudJ::NQD_] = {};
    int NICA = 0, JCOUNT = 0;
    bool LDARK = false;

    auto t_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      // CLOUD_JX mutates CLF in place, so reset it from CLDFRW each call,
      // matching how the normal SZA-scan loop resets it per iteration.
      for (int L = 0; L < LTOP; ++L) {
        CLF[L] = CLDFRW[L];
      }
      rc = 0;
      engine.cloud_jx(
          U0, SZA, RFL, SOLF, LPRTJ,
          PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
          LWP, IWP, REFFL, REFFI, CLF, CLDIW, CLDCOR_in,
          AERSP, NDXAER,
          L1_, CloudJ::AN_, CloudJ::JVN_,
          VALJXX, SKPERD, SWMSQ, OD18,
          IRAN, NICA, JCOUNT, LDARK,
          WTQCA, rc);
      if (rc != 0) {
        std::cerr << "Error: CLOUD_JX failed with rc=" << rc
                   << " at benchmark iteration " << iter << "\n";
        return 1;
      }
    }
    auto t_end = std::chrono::high_resolution_clock::now();

    double total_s =
        std::chrono::duration<double>(t_end - t_start).count();
    double per_call_ms = (total_s / static_cast<double>(bench_iters)) * 1000.0;
    double throughput = static_cast<double>(bench_iters) / total_s;

    std::cout << "BENCH: " << bench_iters << " iters, total="
              << std::fixed << std::setprecision(6) << total_s
              << " s, per-call=" << std::setprecision(6) << per_call_ms
              << " ms, throughput=" << std::setprecision(1) << throughput
              << " calls/s" << std::endl;

    return 0;
  }

  // =====================================================================
  // Print banner (matches Fortran)
  // =====================================================================
  std::cout << ">>>begin Cloud-J v8.0 Standalone" << std::endl;

  // =====================================================================
  // SZA scan loop (matches Fortran: 3 SZAs for test)
  // =====================================================================
  int SZAscan[3] = {0, 30, 60};

  for (int I = 0; I < 3; ++I) {
    int NSZA = SZAscan[I];
    double SZA = static_cast<double>(NSZA);

    // Reset CLF each SZA iteration
    for (int L = 0; L < LTOP; ++L) {
      CLF[L] = CLDFRW[L];
    }

    int IRAN = 1;
    double SOLF = 1.0;
    double U0 = std::cos(SZA * CloudJ::CPI180);

    // ANGLES for OSA: EMU(1:4) + U0
    double ANGLES[5];
    ANGLES[0] = CloudJ::EMU[0];
    ANGLES[1] = CloudJ::EMU[1];
    ANGLES[2] = CloudJ::EMU[2];
    ANGLES[3] = CloudJ::EMU[3];
    ANGLES[4] = U0;

    // Compute OSA and set RFL (surface reflectivity)
    // RFL layout: column-major [5][WW], where RFL[J + 5*K] = albedo for angle J, wvl K
    double RFL[5 * WW] = {};
    const CloudJ::CloudJState& state = engine.get_state();

    for (int K = 0; K < CloudJ::NS2; ++K) {
      double WAVEL = state.WL[K];
      double OSA_dir[5] = {};
      CloudJ::OSA::FJX_OSA(WAVEL, WIND, CHLR, ANGLES, OSA_dir);

      for (int J = 0; J < 5; ++J) {
        RFL[J + 5 * K] = OSA_dir[J];
        // Override OSA with read-in ALBEDO values (same as Fortran test)
        RFL[J + 5 * K] = ALBEDO[J];
      }
    }

    // Set LPRTJ = true to trigger PHOTO_JX printing (matches Fortran)
    bool LPRTJ = true;

    // Output arrays
    double VALJXX[L_ * CloudJ::JVN_] = {};
    double SKPERD[(CloudJ::S_ + 2) * L1_] = {};
    double SWMSQ[6] = {};
    double OD18[L1_] = {};
    double WTQCA[CloudJ::NQD_] = {};
    int NICA = 0, JCOUNT = 0;
    bool LDARK = false;

    // Print SZA info (matches Fortran diagnostic)
    if (LPRTJ) {
      std::cout << std::fixed;
      std::cout << "SZA SOLF U0 albedo"
                << std::setw(8) << std::setprecision(3) << SZA
                << std::setw(8) << std::setprecision(5) << SOLF
                << std::setw(8) << std::setprecision(5) << U0
                << std::setw(8) << std::setprecision(5) << RFL[4 + 5 * (CloudJ::W_ - 1)]
                << std::endl;
      CloudJ::PhotoJX::JP_ATM0(PPP, TTT, DDD, OOO, ZZZ, L_);
      std::cout << " wvl  albedo u1:u4 & u0" << std::endl;
      for (int K = 0; K < CloudJ::NS2; ++K) {
        std::cout << std::setw(5) << (K + 1)
                  << std::fixed << std::setw(8) << std::setprecision(1) << state.WL[K];
        for (int J = 0; J < 5; ++J) {
          std::cout << std::fixed << std::setw(8) << std::setprecision(4) << RFL[J + 5 * K];
        }
        std::cout << std::endl;
      }
    }

    // Call CLOUD_JX
    rc = 0;
    engine.cloud_jx(
        U0, SZA, RFL, SOLF, LPRTJ,
        PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
        LWP, IWP, REFFL, REFFI, CLF, CLDIW, CLDCOR_in,
        AERSP, NDXAER,
        L1_, CloudJ::AN_, CloudJ::JVN_,
        VALJXX, SKPERD, SWMSQ, OD18,
        IRAN, NICA, JCOUNT, LDARK,
        WTQCA, rc);

    if (rc != 0) {
      std::cerr << "Error: CLOUD_JX failed with rc=" << rc << "\n";
      return 1;
    }

    // Print J-values (since PHOTO_JX print not implemented in C++ port)
    // Print J-values (since PHOTO_JX print not implemented in C++ port)
    // Match Fortran format: "Fast-J ----J-values----"
    int NJX = state.NJX;
    std::cout << " Fast-J ----J-values----" << std::endl;
    // Header line with species titles
    std::cout << " L=  ";
    for (int j = 0; j < NJX; ++j) {
      // Fortran uses 72(a6,3x) format
      std::string title = state.TITLEJX[j];
      // Trim trailing spaces
      while (!title.empty() && title.back() == ' ') title.pop_back();
      std::cout << std::setw(6) << std::left << title << "   ";
    }
    std::cout << std::right << std::endl;

    // Print levels L_ down to 1 (Fortran 1-based)
    // VALJXX layout (column-major): VALJXX[l + L_ * j]
    for (int l = L_ - 1; l >= 0; --l) {
      std::cout << std::setw(3) << (l + 1);
      for (int j = 0; j < NJX; ++j) {
        std::cout << format_fortran_e9_2(VALJXX[l + L_ * j]);
      }
      std::cout << std::endl;
    }

  } // end SZA scan

  return 0;
}
