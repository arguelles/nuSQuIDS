#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids {

struct LVParameters {
  gsl_complex c_emu;
  gsl_complex c_mutau;
};

class nuSQUIDSLV: public nuSQUIDS {
  private:
    bool lv_parameters_set = false;
    int n_ = 1;
    LVParameters c_params;
    squids::SU_vector LVP;
    std::vector<squids::SU_vector> LVP_evol;

    void AddToPreDerive(double x);
    void AddToWriteHDF5(hid_t hdf5_loc_id) const;
    void AddToReadHDF5(hid_t hdf5_loc_id);
    squids::SU_vector HI(unsigned int ie, unsigned int irho) const;
    void initialize_LVP_evol_at_each_energy_node();

  public:
    nuSQUIDSLV() {}
    nuSQUIDSLV(marray<double,1> E_vector,unsigned int numneu, NeutrinoType NT = both, bool iinteraction = false, std::shared_ptr<CrossSectionLibrary> ncs = nullptr);
    nuSQUIDSLV(std::string hdf5_filename, std::string grp = "/", std::shared_ptr<InteractionStructure> int_struct = nullptr);
    nuSQUIDSLV(unsigned int numneu, NeutrinoType NT = neutrino);
    
    void Set_LV_OpMatrix(LVParameters &lv_params);
    void Set_LV_OpMatrix(gsl_matrix_complex *cmatrix);
    void Set_LV_Operator(squids::SU_vector op);
    void Set_LV_EnergyPower(int n);
    void Set_MixingAngle(unsigned int i, unsigned int j, double angle);
    void Set_CPPhase(unsigned int i, unsigned int j, double angle);
    void dump_probabilities() const;

    void Set_initial_state(const marray<double, 1> &ini_state, Basis basis = flavor);
    void Set_initial_state(const marray<double, 2> &ini_state, Basis basis = flavor);
    void Set_initial_state(const marray<double, 3> &ini_state, Basis basis = flavor);
};

class nuSQUIDSLVAtm : public nuSQUIDSAtm<nuSQUIDSLV> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSLV>::nuSQUIDSAtm;

    // nuSQUIDSLVAtm has an attribute nusq_array that contains one nuSQUIDSLV object
    // for each zenith. This array can be obtained with this->GetnuSQuIDS(). The Setter
    // methods have to be wrapped such that the Setter of each element of that array is called

    void Set_LV_OpMatrix(LVParameters &lv_params){
      for(nuSQUIDSLV& nsq : this->GetnuSQuIDS()) nsq.Set_LV_OpMatrix(lv_params);
    }

    void Set_LV_OpMatrix(gsl_matrix_complex *cmatrix){
      for(nuSQUIDSLV& nsq : this->GetnuSQuIDS()) nsq.Set_LV_OpMatrix(cmatrix);
    }

    void Set_LV_Operator(squids::SU_vector op){
      for(nuSQUIDSLV& nsq : this->GetnuSQuIDS()) nsq.Set_LV_Operator(op);
    }

    void Set_LV_EnergyPower(int n){
      for(nuSQUIDSLV& nsq : this->GetnuSQuIDS()) nsq.Set_LV_EnergyPower(n);
    }

}; 
} // namespace nusquids

#endif // nusquidslv_H