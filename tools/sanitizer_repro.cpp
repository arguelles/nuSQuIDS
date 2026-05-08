// Minimal CPU-only repro for UBSan. Exercises the nuSQUIDSAtm
// constructor's emplace_back + nuSQUIDS move-ctor path that grew
// the bug surface this PR aims to fix. Built with -fsanitize=undefined.
#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>

int main() {
  squids::Const units;

  auto costh   = nusquids::linspace(-1.0, -0.1, 10);                     // ncz=10 forces extra realloc passes
  auto e_range = nusquids::logspace(1.0e2 * units.GeV, 1.0e6 * units.GeV, 40);

  nusquids::nuSQUIDSAtm<> nsq(costh, e_range, 3, nusquids::both, /*interactions=*/true);

  // Default-constructed members would hit UB if read before write.
  nsq.Set_MixingAngle(0, 1, 0.563942);
  nsq.Set_MixingAngle(0, 2, 0.154085);
  nsq.Set_MixingAngle(1, 2, 0.785398);
  nsq.Set_SquareMassDifference(1, 7.65e-05);
  nsq.Set_SquareMassDifference(2, 0.00247);

  // Flat unit flux for all flavors / both nu types.
  auto ncz = costh.extent(0);
  auto ne  = e_range.extent(0);
  nusquids::marray<double, 4> inistate{ncz, ne, 2u, 3u};
  std::fill(inistate.begin(), inistate.end(), 1.0);
  nsq.Set_initial_state(inistate, nusquids::flavor);

  nsq.Set_ProgressBar(false);   // intentional: probe value-init of progressbar_loop
  nsq.Set_rel_error(1e-6);
  nsq.Set_abs_error(1e-6);

  std::cout << "Evolving..." << std::flush;
  nsq.EvolveState();
  std::cout << " done.\n";

  // Force an EvalFlavor read so any indeterminate field touched late surfaces.
  double v = nsq.EvalFlavor(1, costh[5], e_range[20], 0);
  std::cout << "nu_mu @ midpoint = " << v << "\n";
  return 0;
}
