#ifndef _solar_probabilities_
#define _solar_probabilities_

#include <SQuIDS/const.h>
#include <SQuIDS/SUNalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_eigen.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>

#include "SolarModel.h"

/*
template<typename FunctionType>
double integrate(FunctionType f, double a, double b){
	double (*wrapper)(double,void*)=[](double x, void* params){
		FunctionType& f=*static_cast<FunctionType*>(params);
		return(f(x));
	};

	gsl_integration_workspace* ws=gsl_integration_workspace_alloc(10000);
	double result, error;
	gsl_function F;
	F.function = wrapper;
	F.params = &f;

	gsl_integration_qags(&F, a, b, 0, 1e-2, 10000, ws, &result, &error);
	gsl_integration_workspace_free(ws);

	return(result);
}
*/

class SOP: public nusquids::nuSQUIDS {
  private:
    /// \brief solar model object
    std::shared_ptr<SolarModel> solar_model = nullptr;
//    /// \brief number of neutrino flavors.
//    const unsigned int numneu = 3;
    /// \brief flavors
    enum flavors {nue = 0,numu = 1,nutau = 2};
    enum specifies {neutrino = 0,antineutrino = 1};
  protected:
    squids::SU_vector Hamiltonian(double E, double r) const;
  public:
    SOP(unsigned int numneu): nusquids::nuSQUIDS(numneu){}
    SOP(): SOP(3){}

    void SetSolarModel(std::shared_ptr<SolarModel> solar_model_){
      solar_model = solar_model_;
    }

    std::shared_ptr<SolarModel> GetSolarModel(){
      if (solar_model == nullptr)
        throw std::runtime_error("No solar model set.");
      return solar_model;
    }

    double SolarOscillationProbability(double E,double r) const;
    double PeeSquare(double E) const;
    double RadialIntegratedFluxes(double E) const;
};

#endif
