/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS

See https://arxiv.org/abs/2007.00068 for model description

Tom Stuttard, Mikkel Jensen, Martin Langgard Ravn (Niels Bohr Institute)
*/


#ifndef NUSQUIDSDECOH_H
#define NUSQUIDSDECOH_H

#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids{


/*
  nuSQuIDS extended to include decoherence
*/

class nuSQUIDSDecoh: public nuSQUIDS {
  public:
    /// \brief Neutrino decoherence models.
    enum DecoherenceModel {RandomizePhase = 0, RandomizeState = 1, NeutrinoLoss = 2};
    //Default void constructor (never used but required by boost python, specifically RegisterBasicNuSQuIDSPythonBindings, not sure why).
    nuSQUIDSDecoh() {}

    //TODO move constructor

    /// Multiple energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSDecoh(marray<double,1> E_vector,unsigned int numneu, NeutrinoType NT = both, bool iinteraction = false, std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
      : nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs)
    { init(); }

    /// Single energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSDecoh(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu, NT)
    { init(); }

    // Load from serialized file (not yet implemented, but required to exist for pybindings)
    nuSQUIDSDecoh(std::string hdf5_filename, std::string grp = "/", std::shared_ptr<InteractionStructure> int_struct = nullptr)
    {
      throw std::runtime_error("nuSQUIDSDecoh::Error::HDF5 serialization not implemented");
    }
  
    // Enable/disable the decoherence D[rho] operator
    void EnableDecoherence(bool enable);

    // Set the Gamma decoherence matrix elements individually
    // `standard_gell_mann` indicates the choice of basis vectors in which the matrix defined
    void Set_DecoherenceGammaMatrix(const marray<double,2>& dmat, bool standard_gell_mann=true);

    // Set the diagonal elements of the Gamma matrix (many models only consider these)
    void Set_DecoherenceGammaMatrixDiagonal(const marray<double,1>& dmat, bool standard_gell_mann=true);

    // Set decoherence models according to specific pre-build models
    void Set_DecoherenceGammaMatrix(DecoherenceModel decoh_model, double gamma);

    // Get the current value of the decoherence Gamma matrix
    marray<double,2> Get_DecoherenceGammaMatrix() const;

    // Get/set the energy dependence using the form:  Gamma = Gamma_0 * (E/E_0) ^ n_energy
    // See https://arxiv.org/abs/1803.04438.pdf equation 2.13
    void Set_DecoherenceGammaEnergyDependence(double n_energy);
    double Get_DecoherenceGammaEnergyDependence() const;

    // Get/set energy scale (E0) relative to which decoherence Gamma matrix is defined
    void Set_DecoherenceGammaEnergyScale(double energy);
    double Get_DecoherenceGammaEnergyScale() const;


  protected:

    // Common initialisation tasks
    void init();

    // Handle mapping based on generator choice
    unsigned int MapBasisVectorConventions(unsigned int i);


  private:

    // Store the number of basis vectors
    unsigned short num_basis_vectors;

    // The matrix containing the Gamma coefficients in D[rho]
    // This is defined in the SU(N) basis
    //TODO ref my paper
    std::unique_ptr<gsl_matrix_complex> decoherence_gamma_matrix;

    // Index of Gamma energy dependence
    double gamma_energy_dep_index;

    // Energy scale
    double gamma_energy_scale; // E0

    // Function to return D[rho] (result of decohernce operator on rho)
    squids::SU_vector DRho(unsigned int ei,unsigned int index_rho) const;

    // Overrides nuSQuIDS Interactions Rho to add DRho terms 
    squids::SU_vector InteractionsRho(unsigned int ei, unsigned int irho) const override;

    // Units container
    squids::Const units;

    // bool to enable decoherence
    bool enable_decoherence;
};


/*
  nuSQUIDSAtm extended to include decoherence
*/

class nuSQUIDSDecohAtm : public nuSQUIDSAtm<nuSQUIDSDecoh> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSDecoh>::nuSQUIDSAtm;

    // TODO Copy and move constructors...

    // Wrap all the getters/setters
    void Set_Debug(bool debug) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.Set_Debug(debug);
    } 

    void EnableDecoherence(bool enable) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.EnableDecoherence(enable);
    }

    void Set_DecoherenceGammaMatrix(const marray<double,2>& dmat) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.Set_DecoherenceGammaMatrix(dmat);
    }

    void Set_DecoherenceGammaMatrixDiagonal(const marray<double,1>& dmat, bool standard_gell_mann=true) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.Set_DecoherenceGammaMatrixDiagonal(dmat, standard_gell_mann);
    }

    marray<double,2> Get_DecoherenceGammaMatrix() {
      return this->GetnuSQuIDS(0).Get_DecoherenceGammaMatrix();
    }

    void Set_DecoherenceGammaEnergyDependence(double n_energy) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.Set_DecoherenceGammaEnergyDependence(n_energy);
    }

    double Get_DecoherenceGammaEnergyDependence() {
      return this->GetnuSQuIDS(0).Get_DecoherenceGammaEnergyDependence();
    }

    void Set_DecoherenceGammaEnergyScale(double energy) {
      for(nuSQUIDSDecoh& nsq : this->GetnuSQuIDS()) nsq.Set_DecoherenceGammaEnergyScale(energy);
    }

    double Get_DecoherenceGammaEnergyScale() {
      return this->GetnuSQuIDS(0).Get_DecoherenceGammaEnergyScale();
    }

};



} // close namespace

#endif //NUSQUIDSDECOH_H
