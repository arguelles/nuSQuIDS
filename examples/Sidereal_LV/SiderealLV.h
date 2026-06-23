#ifndef NUSQUIDS_SIDEREAL_LV_H
#define NUSQUIDS_SIDEREAL_LV_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <nuSQuIDS/nuSQuIDS.h>

// Implementation of neutrino oscillations with Sidereal Lorentz Violation (LV)
// from the Standard Model Extension (SME), minimal sector.
// Follows https://doi.org/10.1103/PhysRevD.85.096005.
//
// SME coefficients:
//   fa[i][j][k]       - aT parameters (mass dim 3, CPT-odd), units of eV
//   fc[i][j][k1][k2]  - cT parameters (mass dim 4, CPT-even), dimensionless
//
// Flavour convention: 0=e, 1=mu, 2=tau
// Spatial coord convention: 0=X, 1=Y, 2=Z (SME Sun-centred frame)

namespace nusquids {

class nuSQUIDSSiderealLV : public nuSQUIDS {
 private:
  // SME coefficients (upper triangle in flavour space, i <= j)
  double fa[3][3][3];       // aT [flvi][flvj][coord], in eV
  double fc[3][3][3][3];    // cT [flvi][flvj][coord1][coord2], dimensionless

  // Neutrino direction factors (NX, NY, NZ) in the SME frame
  double fN[3];

  // Sidereal time and detector geometry
  double fTime;  // local sidereal time, in hours
  double fChi;   // detector colatitude (90 - latitude), in degrees

  // Per-energy LV vectors in mass basis (interaction picture)
  // aT: CPT-odd (sign flips for antineutrinos)
  // cT: CPT-even (same sign for neutrinos and antineutrinos)
  std::vector<squids::SU_vector> LV_aT_base;
  std::vector<squids::SU_vector> LV_cT_base;
  std::vector<squids::SU_vector> LV_aT_evol;
  std::vector<squids::SU_vector> LV_cT_evol;

  bool dirty; // triggers recomputation of LV_base vectors

  static constexpr double kSiderealDayHours = 23.9344696;
  static constexpr double kOmegaSidereal    = 2.0 * M_PI / kSiderealDayHours;

  // Build LV base vectors (mass basis, not yet evolved) for all energy bins.
  // Called once whenever SME parameters, fTime, or fN change.
  void ComputeLVBase() {
    double st  = std::sin(kOmegaSidereal * fTime);
    double ct  = std::cos(kOmegaSidereal * fTime);
    double s2t = 2.0 * st * ct;
    double c2t = ct * ct - st * st;

    for (unsigned int ei = 0; ei < ne; ei++) {
      double E = E_range[ei]; // energy in eV (SQuIDS internal units)

      gsl_matrix_complex* M_aT = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* M_cT = gsl_matrix_complex_calloc(3, 3);

      for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
          // aT contributions (sign = +1 for neutrinos; will be negated for antineutrinos)
          double C0  = -fa[i][j][2] * fN[2];
          double As0 =  fa[i][j][0] * fN[1] - fa[i][j][1] * fN[0];
          double Ac0 = -(fa[i][j][0] * fN[0] + fa[i][j][1] * fN[1]);
          double val_aT = C0 + As0 * st + Ac0 * ct;

          // cT contributions (CPT-even, no sign flip)
          double As1 =  2.0*fN[1]*fN[2]*fc[i][j][0][2]
                      - 2.0*fN[0]*fN[2]*fc[i][j][1][2];
          double Ac1 = -2.0*fN[0]*fN[2]*fc[i][j][0][2]
                      - 2.0*fN[1]*fN[2]*fc[i][j][1][2];
          double Bs1 =  fN[0]*fN[1]*(fc[i][j][0][0] - fc[i][j][1][1])
                      - (fN[0]*fN[0] - fN[1]*fN[1])*fc[i][j][0][1];
          double Bc1 = -0.5*(fN[0]*fN[0] - fN[1]*fN[1])
                           *(fc[i][j][0][0] - fc[i][j][1][1])
                      - 2.0*fN[0]*fN[1]*fc[i][j][0][1];
          double val_cT = E * (As1*st + Ac1*ct + Bs1*s2t + Bc1*c2t);

          gsl_complex z_aT = gsl_complex_rect(val_aT, 0.0);
          gsl_complex z_cT = gsl_complex_rect(val_cT, 0.0);

          gsl_matrix_complex_set(M_aT, i, j, z_aT);
          gsl_matrix_complex_set(M_cT, i, j, z_cT);
          if (i != j) {
            gsl_matrix_complex_set(M_aT, j, i, z_aT);
            gsl_matrix_complex_set(M_cT, j, i, z_cT);
          }
        }
      }

      LV_aT_base[ei] = squids::SU_vector(M_aT);
      LV_cT_base[ei] = squids::SU_vector(M_cT);
      LV_aT_base[ei].RotateToB1(params);
      LV_cT_base[ei].RotateToB1(params);

      gsl_matrix_complex_free(M_aT);
      gsl_matrix_complex_free(M_cT);
    }
    dirty = false;
  }

  void AddToPreDerive(double x) {
    if (dirty) ComputeLVBase();
    double dt = x - Get_t_initial();
    for (unsigned int ei = 0; ei < ne; ei++) {
      LV_aT_evol[ei] = LV_aT_base[ei].Evolve(H0_array[ei], dt);
      LV_cT_evol[ei] = LV_cT_base[ei].Evolve(H0_array[ei], dt);
    }
  }

  squids::SU_vector HI(unsigned int ei, unsigned int irho) const {
    squids::SU_vector potential = nuSQUIDS::HI(ei, irho);
    bool is_antineu = (irho == 1 && NT == both) || NT == antineutrino;
    if (is_antineu)
      potential += (-1.0) * LV_aT_evol[ei] + LV_cT_evol[ei];
    else
      potential += LV_aT_evol[ei] + LV_cT_evol[ei];
    return potential;
  }

  // Validate flavour pair (i <= j required). Returns false and warns on error.
  static bool ValidateAndOrder(int& flvi, int& flvj, const char* name) {
    if (flvi > flvj) {
      std::cerr << "WARNING: " << name << " flavour indices swapped ("
                << flvj << "," << flvi << ")." << std::endl;
      std::swap(flvi, flvj);
    }
    if (flvi < 0 || flvj > 2) {
      std::cerr << "WARNING: " << name << "_" << flvi << flvj
                << " invalid for 3 neutrinos. Doing nothing." << std::endl;
      return false;
    }
    return true;
  }

  // Validate spatial index in [0,2].
  static bool ValidateCoord(int coord, const char* coordName, const char* name) {
    if (coord < 0 || coord > 2) {
      std::cerr << "WARNING: " << coordName << "=" << coord
                << " out of range [0,2] for " << name << ". Doing nothing." << std::endl;
      return false;
    }
    return true;
  }

 public:
  nuSQUIDSSiderealLV(marray<double,1> E_range_, unsigned int numneu_,
                     NeutrinoType NT_, bool interaction_)
      : nuSQUIDS(E_range_, numneu_, NT_, interaction_),
        fa{}, fc{}, fN{0.0, 0.0, 0.0}, fTime(0.0), fChi(0.0), dirty(true)
  {
    if (numneu_ != 3)
      throw std::runtime_error("nuSQUIDSSiderealLV: only 3 flavours supported");
    LV_aT_base.resize(ne, squids::SU_vector(nsun));
    LV_cT_base.resize(ne, squids::SU_vector(nsun));
    LV_aT_evol.resize(ne, squids::SU_vector(nsun));
    LV_cT_evol.resize(ne, squids::SU_vector(nsun));
  }

  // ---- SME parameter setters ----

  void SetA(int flvi, int flvj, int coord, double val) {
    if (!ValidateAndOrder(flvi, flvj, "aT")) return;
    if (!ValidateCoord(coord, "coord", "aT")) return;
    fa[flvi][flvj][coord] = val;
    dirty = true;
  }

  void SetC(int flvi, int flvj, int coord1, int coord2, double val) {
    if (!ValidateAndOrder(flvi, flvj, "cT")) return;
    if (!ValidateCoord(coord1, "coord1", "cT")) return;
    if (!ValidateCoord(coord2, "coord2", "cT")) return;
    fc[flvi][flvj][coord1][coord2] = val;
    dirty = true;
  }

  double GetA(int flvi, int flvj, int coord) {
    if (!ValidateAndOrder(flvi, flvj, "aT")) return 0.0;
    if (!ValidateCoord(coord, "coord", "aT")) return 0.0;
    return fa[flvi][flvj][coord];
  }

  double GetC(int flvi, int flvj, int coord1, int coord2) {
    if (!ValidateAndOrder(flvi, flvj, "cT")) return 0.0;
    if (!ValidateCoord(coord1, "coord1", "cT")) return 0.0;
    if (!ValidateCoord(coord2, "coord2", "cT")) return 0.0;
    return fc[flvi][flvj][coord1][coord2];
  }

  // ---- Geometry and time setters ----

  // Set the neutrino direction from zenith and azimuth (degrees).
  // Requires SetColatitude() to be called first.
  void SetNeutrinoDirection(double zenith_deg, double azimuth_deg) {
    double chi = fChi        * M_PI / 180.0;
    double zen = zenith_deg  * M_PI / 180.0;
    double azi = azimuth_deg * M_PI / 180.0;
    fN[0] =  std::cos(chi)*std::sin(zen)*std::cos(azi) + std::sin(chi)*std::cos(zen);
    fN[1] =  std::sin(zen)*std::sin(azi);
    fN[2] = -std::sin(chi)*std::sin(zen)*std::cos(azi) + std::cos(chi)*std::cos(zen);
    dirty = true;
  }

  // Set local sidereal time in hours.
  void SetTimeHours(double hours) {
    fTime = hours;
    dirty = true;
  }

  // Set detector colatitude (chi = 90 - latitude) in degrees.
  void SetColatitude(double chi_deg) {
    fChi = chi_deg;
  }

  // Set detector colatitude from geographic coordinates (positive = North).
  void SetColatitude(double deg, double min, double sec = 0.0) {
    double latitude = deg + min / 60.0 + sec / 3600.0;
    fChi = 90.0 - latitude;
  }

  double GetColatitude() const { return fChi; }

  // Override mixing-angle/phase setters to invalidate LV rotation cache.
  void Set_MixingAngle(unsigned int i, unsigned int j, double angle) {
    nuSQUIDS::Set_MixingAngle(i, j, angle);
    dirty = true;
  }

  void Set_CPPhase(unsigned int i, unsigned int j, double phase) {
    nuSQUIDS::Set_CPPhase(i, j, phase);
    dirty = true;
  }

  // Preserve dirty flag across Set_initial_state calls (nuSQuIDS resets params).
  void Set_initial_state(const marray<double,1>& v, Basis b) {
    bool d = dirty;
    nuSQUIDS::Set_initial_state(v, b);
    dirty = d;
  }
  void Set_initial_state(const marray<double,2>& v, Basis b) {
    bool d = dirty;
    nuSQUIDS::Set_initial_state(v, b);
    dirty = d;
  }
  void Set_initial_state(const marray<double,3>& v, Basis b) {
    bool d = dirty;
    nuSQUIDS::Set_initial_state(v, b);
    dirty = d;
  }
};

} // namespace nusquids

#endif // NUSQUIDS_SIDEREAL_LV_H
