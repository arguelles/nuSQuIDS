// Example: P(numu->numu) with Sidereal LV vs Right Ascension for IceCube.
//
// Scans Right Ascension (0-360 deg), averaging over Declination,
// using IceCube geometry and Earth (PREM) matter potential.
// Single fixed energy node (no energy loop, no spline interpolation).
//
// SME parameter: aT[mu][tau][X] = 2e-23 eV (only non-zero LV coefficient)
// Energy: 1 TeV
// Detector: IceCube [89deg 59' 24'' S]
//
// Compile:  make examples/Sidereal_LV/sidereal_lv
// Output:   sidereal_lv_flux.txt  (RA[deg]  P_std  P_lv)

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

#include <nuSQuIDS/nuSQuIDS.h>
#include "SiderealLV.h"

using namespace nusquids;

int main() {
  squids::Const units;

  // ---- IceCube geometry ----
  // IceCube is at the South Pole: 89deg 59' 24'' S
  // colatitude = 90 + (deg + min/60 + sec/3600) for Southern latitudes
  const double ic_deg     = 89.0, ic_min = 59.0, ic_sec = 24.0;
  const double chi_ic_deg = 90.0 + (ic_deg + ic_min/60.0 + ic_sec/3600.0);
  const double lat_ic     = (90.0 - chi_ic_deg) * M_PI / 180.0; // latitude in rad (~-pi/2)

  // ---- Physics parameters ----
  const double E       = 1.0 * units.TeV;    // fixed energy, single node
  const double aT_val  = 2.0e-23 * units.GeV; // aT[mu][tau][X], only non-zero LV term (GeV scale, as in OscProb convention)
  const double dm21    = 7.5e-5  * units.eV * units.eV;
  const double dm31    = 2.457e-3 * units.eV * units.eV;
  const double R_Earth = 6371.0 * units.km;

  // ---- Scan grid ----
  const int N_RA  = 360;  // Right Ascension bins (1 deg steps)
  const int N_DEC = 90;   // Declination bins (2 deg steps)

  // Single energy node: no energy loop, no spline interpolation.
  marray<double,1> Erange({1});
  Erange[0] = E;

  marray<double,2> ini({1, 3});
  ini[0][0] = 0.0; ini[0][1] = 1.0; ini[0][2] = 0.0; // pure numu

  // ---- Build objects ONCE and reuse them across the scan ----
  auto earth = std::make_shared<Earth>(); // loads PREM profile once

  nuSQUIDS nus_std(Erange, 3, neutrino, false);
  nus_std.Set_MixingAngle(0, 1, 0.5836);
  nus_std.Set_MixingAngle(0, 2, 0.1496);
  nus_std.Set_MixingAngle(1, 2, 0.8587);
  nus_std.Set_CPPhase(0, 2, -1.601);
  nus_std.Set_SquareMassDifference(1, dm21);
  nus_std.Set_SquareMassDifference(2, dm31);
  nus_std.Set_Body(earth);
  nus_std.Set_rel_error(1.0e-6);
  nus_std.Set_abs_error(1.0e-6);
  nus_std.Set_ProgressBar(false);

  nuSQUIDSSiderealLV nus_lv(Erange, 3, neutrino, false);
  nus_lv.Set_MixingAngle(0, 1, 0.5836);
  nus_lv.Set_MixingAngle(0, 2, 0.1496);
  nus_lv.Set_MixingAngle(1, 2, 0.8587);
  nus_lv.Set_CPPhase(0, 2, -1.601);
  nus_lv.Set_SquareMassDifference(1, dm21);
  nus_lv.Set_SquareMassDifference(2, dm31);
  nus_lv.SetA(1, 2, 0, aT_val);      // aT[mu][tau][X] = 2e-23 eV (only non-zero term)
  nus_lv.SetColatitude(chi_ic_deg);  // IceCube South Pole colatitude
  nus_lv.SetTimeHours(0.0);
  nus_lv.Set_Body(earth);
  nus_lv.Set_rel_error(1.0e-6);
  nus_lv.Set_abs_error(1.0e-6);
  nus_lv.Set_ProgressBar(false);

  std::ofstream file("sidereal_lv_flux.txt");
  file << "# RA[deg]  P_std(numu->numu)  P_lv(numu->numu)\n";

  for (int i_ra = 0; i_ra < N_RA; i_ra++) {
    double alpha = (i_ra + 0.5) * 360.0 / N_RA;  // RA in deg, bin centre

    double sum_std = 0.0;
    double sum_lv  = 0.0;
    int    count   = 0;

    for (int i_dec = 0; i_dec < N_DEC; i_dec++) {
      double delta_deg = -90.0 + (i_dec + 0.5) * 180.0 / N_DEC;  // DEC bin centre
      double delta     = delta_deg * M_PI / 180.0;

      // Zenith angle at IceCube for this (RA, DEC) — same formula as SiderealLIV2D.C
      double azimuth = M_PI + alpha * M_PI / 180.0;
      double cos_zen = std::sin(lat_ic) * std::sin(delta)
                     + std::cos(lat_ic) * std::cos(delta) * std::cos(azimuth);
      cos_zen = std::max(-1.0, std::min(1.0, cos_zen));
      double zenith_deg = std::acos(cos_zen) * 180.0 / M_PI;

      // Only upgoing neutrinos cross the Earth (zenith > 90 deg)
      if (cos_zen >= 0.0) continue;

      // Chord length through Earth
      double L = -2.0 * R_Earth * cos_zen;  // cos_zen < 0, so L > 0
      auto track = std::make_shared<Earth::Track>(L);

      // ---- Standard nuSQuIDS (no LV) ----
      nus_std.Set_Track(track);
      nus_std.Set_initial_state(ini, flavor);
      nus_std.EvolveState();
      sum_std += nus_std.EvalFlavorAtNode(1, 0, 0); // numu, energy node 0, neutrino

      // ---- Sidereal LV nuSQuIDS ----
      nus_lv.SetNeutrinoDirection(zenith_deg, alpha + 180.0); // azimuth = pi + alpha
      nus_lv.Set_Track(track);
      nus_lv.Set_initial_state(ini, flavor);
      nus_lv.EvolveState();
      sum_lv += nus_lv.EvalFlavorAtNode(1, 0, 0);

      count++;
    }

    if (count > 0)
      file << alpha << "  " << sum_std / count << "  " << sum_lv / count << "\n";

    if (i_ra % 10 == 0)
      std::cout << "RA = " << alpha << " deg  (" << i_ra+1 << "/" << N_RA << ")\r" << std::flush;
  }

  std::cout << "\nDone! Output written to sidereal_lv_flux.txt" << std::endl;
  return 0;
}
