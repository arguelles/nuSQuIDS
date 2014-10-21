#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQUIDS.h>

using namespace nusquids;

void exercise_se_mode(unsigned int numneu,string NT, std::shared_ptr<Body> body, std::shared_ptr<Track> track){
  nuSQUIDS nus(numneu,NT);
  nus.Set_Track(track);
  nus.Set_Body(body);

  switch (numneu){
    case 3:
      nus.Set(TH12,0.583996);
      nus.Set(TH13,0.148190);
      nus.Set(TH23,0.737324);
      nus.Set(DM21SQ,7.5e-05);
      nus.Set(DM31SQ,0.00257);
      nus.Set(DELTA1,1.);
      break;
    case 4:
      // random values for non standart parameters
      nus.Set(TH12,0.583996);
      nus.Set(TH13,0.148190);
      nus.Set(TH23,0.737324);
      nus.Set(TH14,0.1245);
      nus.Set(TH24,0.5454);
      nus.Set(TH34,0.32974);
      nus.Set(DM21SQ,7.5e-05);
      nus.Set(DM31SQ,0.00257);
      nus.Set(DM41SQ,1.9234);
      nus.Set(DELTA1,1.);
      nus.Set(DELTA2,0.135);
      break;
  }

  // energies
  std::vector<double> test_energies {1.0e-6,1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,1.0e0,1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,1.0e6};

  std::cout << std::setprecision(3);
  std::cout << std::scientific;
  for(int flv = 0; flv < numneu; flv++){
    std::vector<double> ini_state(numneu);
    for (int iflv = 0; iflv < numneu; iflv++)
      ini_state[iflv] = ( iflv==flv ? 1.0 : 0.0 );
      for(double Enu : test_energies){
        nus.Set_initial_state(ini_state,"flavor");
        nus.Set_E(Enu*nus.units.GeV);
        nus.EvolveState();
        std::cout << body->name << " " << flv << " [flv] " << Enu << " [GeV] ";
        for (int i = 0; i < numneu; i++){
          double p = nus.EvalFlavor(i);
          if ( p < 1.0e-8)
            std::cout << 0.0 << " ";
          else
            std::cout << p << " ";
        }
        std::cout << std::endl;
      }
  }
}

int main(){
  // this test checks the neutrino oscillation probability
  // for 3 and 4 neutrino flavors for neutrinos/antineutrinos/both
  // and the single and multiple energy modes

  // units
  Const units;

  /*
  // vacuum
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(0.0,100.0*units.km);

  exercise_se_mode(3,"neutrino",vacuum,track_vac);
  exercise_se_mode(3,"antineutrino",vacuum,track_vac);
  exercise_se_mode(4,"neutrino",vacuum,track_vac);
  exercise_se_mode(4,"antineutrino",vacuum,track_vac);
  */
  // constant density
  double density = 20.0; // gr/cm^3
  double ye = 0.3;
  std::shared_ptr<ConstantDensity> constdens = std::make_shared<ConstantDensity>(density,ye);
  std::shared_ptr<ConstantDensity::Track> track_constdens = std::make_shared<ConstantDensity::Track>(0.0,700.0*units.km);

  exercise_se_mode(3,"neutrino",constdens,track_constdens);
  exercise_se_mode(3,"antineutrino",constdens,track_constdens);
  exercise_se_mode(4,"neutrino",constdens,track_constdens);
  exercise_se_mode(4,"antineutrino",constdens,track_constdens);

  /*
  // earth
  std::shared_ptr<Earth> earth = std::make_shared<Earth>();
  std::shared_ptr<Earth::Track> earth_track = std::make_shared<Earth::Track>(0.0,4000.0*units.km,4000.0*units.km);

  exercise_se_mode(3,"neutrino",earth,earth_track);
  exercise_se_mode(3,"antineutrino",earth,earth_track);
  exercise_se_mode(4,"neutrino",earth,earth_track);
  exercise_se_mode(4,"antineutrino",earth,earth_track);

  // earth atmospheric
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> earth_atm_track = std::make_shared<EarthAtm::Track>(acos(-1));

  exercise_se_mode(3,"neutrino",earth_atm,earth_atm_track);
  exercise_se_mode(3,"antineutrino",earth_atm,earth_atm_track);
  exercise_se_mode(4,"neutrino",earth_atm,earth_atm_track);
  exercise_se_mode(4,"antineutrino",earth_atm,earth_atm_track);

  // sun
  std::shared_ptr<Sun> sun = std::make_shared<Sun>();
  std::shared_ptr<Sun::Track> track_sun = std::make_shared<Sun::Track>(0.0,sun->radius*units.km);

  exercise_se_mode(3,"neutrino",sun,track_sun);
  exercise_se_mode(3,"antineutrino",sun,track_sun);
  exercise_se_mode(4,"neutrino",sun,track_sun);
  exercise_se_mode(4,"antineutrino",sun,track_sun);
  */
}
