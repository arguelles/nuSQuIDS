#include <math.h>
#include "nuSQUIDSDecoh.h"

namespace nusquids{

void nuSQUIDSDecoh::EnableDecoherence(bool enable) {
  enable_decoherence = enable;
}


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // Performing calculation in the interaction basi (which is the default)
  // This is the only basis for which splined interpolation between nodes is implemented
  // Could optionally try solving in e.g. the mass basis (Set_Basis(mass)), which would be simpler for debugging, but have discovered an issue with mass/falvor bassi calculation that needs resolving first

  // Currently supports 3 flavors
  if (numneu != 3) {
    std::cout << "WARNING : Code has only been tested for 3 neutrino flavors" << std::endl; 
  }

  // Set the number of basis vectors in SU(N)
  // This is N^2, where N is the num neutrino flavors (generators + identity matrix)
  num_basis_vectors = pow(numneu, 2); //TODO already a nuSQuiDS variable somewhere specifiying this?

  // Enable decoherence term in the SQuIDS numerical solver
  EnableDecoherence(true);

  // Initialise the decoherence matrix 
  // Allocate the memory, and set all values to 0 (e.g. no decoherence)
  // decoherence_gamma_matrix.SetAllComponents(0.);
  decoherence_gamma_matrix.reset( gsl_matrix_complex_alloc(num_basis_vectors, num_basis_vectors) );
  gsl_matrix_complex_set_zero(decoherence_gamma_matrix.get());

  // Choose some defaults
  // Set_DecoherenceBasis(mass);
  Set_DecoherenceGammaEnergyDependence(0.); // Energy-independent
  Set_DecoherenceGammaEnergyScale(1.*units.GeV); // E0 = 1 GeV
  
}


unsigned int nuSQUIDSDecoh::MapBasisVectorConventions(unsigned int i) {
  /*
    SQuIDS decomposes the various operators (such as density matrix, Hamiltonian, etc) into 
    coefficients of the SU(N) basis vectors (see https://arxiv.org/pdf/1412.3832.pdf eqns 15-17), 
    and uses these for the time evolution. This is often called the "SU(N) basis".

    Note that this happens to match well with the Lindblad master equation which is normally 
    expressed in the SU(N) basis. 
    
    However, this conversion depends on your convention for the basis vectors (e.g. identify 
    matrix + generators), and SQuIDS uses general definition that supports SU(N) (e.g. generalised 
    N, not 2 or 3 specifically) but happens to produce basis vectors in a different order (and 
    occaisionally but unimportantly differnt signs) to the "standard" Gell-Mann matrices definition 
    for the SU(3) case (as shown on https://en.wikipedia.org/wiki/Gell-Mann_matrices).

    Therefore, if I specify some values directly in the SU(3) basis (as people typically do with the 
    8x8 decoherence/dissipation operator coefficients), rather than defining a NxN matrix and using 
    SQuIDS's own conversion, then I need to correct the mapping if the user is assuming the "standard" 
    Gell-Mann matrices definition.
  */

  // Checks
  if(i > num_basis_vectors)
      throw std::runtime_error("nuSQUIDSDecoh::Error:Invalid index " + std::to_string(i) + " for SU(" + std::to_string(numneu) + ")" );

  // SU(2)
  // Mapping is unchanged (e.g. SQuIDS version == Dirac matrices)
  if(numneu == 2) { 
    if(i == 0)      return 0;
    else if(i == 1) return 1;
    else if(i == 2) return 2;
    else if(i == 3) return 3;
  }

  // SU(3)
  // Need to map from Gell-Mann to SQuIDS
  else if(numneu == 3) {
    if(i == 0)      return 0;
    else if(i == 1) return 1;
    else if(i == 2) return 3; // -ve
    else if(i == 3) return 4;
    else if(i == 4) return 2;
    else if(i == 5) return 6; // -ve
    else if(i == 6) return 5;
    else if(i == 7) return 7; // -ve
    else if(i == 8) return 8;

  }

  // Error handling
  throw std::runtime_error("nuSQUIDSDecoh::Error:Generator mapping not defined for SU(" + std::to_string(numneu) + ")");
}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(const marray<double,2>& dmat, bool standard_gell_mann) {

  // Check dimensions
  // Should be N^2 x N^2
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
      if( dmat.extent(dim) != num_basis_vectors)
        throw std::runtime_error("nuSQUIDSDecoh::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(num_basis_vectors) + ")" );
  }

  // Loop over matrix elements (2D)
  for( unsigned int i = 0 ; i < num_basis_vectors ; ++i ) {
      for( unsigned int j = 0 ; j < num_basis_vectors ; ++j ) {

          // Create the complex number to write
          gsl_complex c{{ dmat[i][j] , 0.0 }}; //Only using real part right now

          // Handle mapping based on choice of SU(3) generators
          unsigned int m = i;
          unsigned int n = j;
          if(standard_gell_mann) { //TODO handle SU( N != 3 )
              m = this->MapBasisVectorConventions(i);
              n = this->MapBasisVectorConventions(j);
          }

          // Write to the matrix element
          gsl_matrix_complex_set(decoherence_gamma_matrix.get(), m, n, c);

      }
  }

}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrixDiagonal(const marray<double,1>& dmat, bool standard_gell_mann) {

    // Check dimensions
    // Should be N^2 elements
    if( dmat.extent(0) != num_basis_vectors)
        throw std::runtime_error("nuSQUIDSDecoh::Error:Input diagonal decoherence matrix has wrong length (found " + std::to_string(dmat.extent(0))+ ", should be " + std::to_string(num_basis_vectors) + ")" );

    // Create the 2D matrix using the diagonals provided (off-diagonals are zero)
    //TODO Is there a more efficient way to do this? This is called at each step in a fit so is called multiple times.
    marray<double,2> dmat_2d{num_basis_vectors, num_basis_vectors};
    for( unsigned int i = 0 ; i < num_basis_vectors ; ++i ) {
        for( unsigned int j = 0 ; j < num_basis_vectors ; ++j ) {
            dmat_2d[i][j] =  i == j ? dmat[i] : 0. ;
        }
    }

    Set_DecoherenceGammaMatrix(dmat_2d, standard_gell_mann);

}

void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(DecoherenceModel decoh_model, double gamma){
	marray<double,1> dmat{num_basis_vectors};
	for(int i = 0; i < dmat.extent(0); i++)
		dmat[i] = gamma;
	switch(decoh_model){
		case RandomizePhase: dmat[0] = 0, dmat[3] = 0, dmat[8] = 0;
		case RandomizeState: dmat[0] = 0;
		case NeutrinoLoss: ;
	}
	Set_DecoherenceGammaMatrixDiagonal(dmat,false);
}

marray<double,2> nuSQUIDSDecoh::Get_DecoherenceGammaMatrix() const {
    /*
    Get the current value of the decoherence gamma matrix
    Return the value as the marray type, which can be converted to a numpy array in pyhton
    This is not used as part of the main solver, just for user scripts
    */

    marray<double,2> dmat{num_basis_vectors, num_basis_vectors};
    for( unsigned int i = 0 ; i < num_basis_vectors ; ++i ) {
        for( unsigned int j = 0 ; j < num_basis_vectors ; ++j ) {
          dmat[i][j] = GSL_REAL(gsl_matrix_complex_get(decoherence_gamma_matrix.get(),i,j)); //TODO imaginary components
        }
    }

    return dmat;
}


void nuSQUIDSDecoh::Set_DecoherenceGammaEnergyDependence(double n)  {
      gamma_energy_dep_index = n; 
}

double nuSQUIDSDecoh::Get_DecoherenceGammaEnergyDependence() const {
    return gamma_energy_dep_index; 
}

void nuSQUIDSDecoh::Set_DecoherenceGammaEnergyScale(double energy)  {
    gamma_energy_scale = energy; 
}

double nuSQUIDSDecoh::Get_DecoherenceGammaEnergyScale() const {
    return gamma_energy_scale; 
}


squids::SU_vector nuSQUIDSDecoh::DRho(unsigned int ei,unsigned int index_rho) const {

    // Scale the gamma matrix according to the energy dependence
    double energy_dependence = pow( (E_range[ei] / gamma_energy_scale) , gamma_energy_dep_index);

    // Compute the decoherence operator, D[rho]
    auto rho_components = estate[ei].rho[index_rho].GetComponents();
    std::vector<double> DRho_components(num_basis_vectors, 0.);
    for( unsigned int i = 0 ; i < num_basis_vectors ; ++i ) {
        DRho_components[i] = 0.;
        for( unsigned int j = 0 ; j < num_basis_vectors ; ++j ) {
            DRho_components[i] += GSL_REAL(gsl_matrix_complex_get(decoherence_gamma_matrix.get(),i,j)) * rho_components[j] * energy_dependence; //TODO imag
        }
    }
    squids::SU_vector DRho_value(DRho_components);

    return DRho_value;

}


squids::SU_vector nuSQUIDSDecoh::InteractionsRho(unsigned int ei, unsigned int irho) const {
  if (enable_decoherence)
    return nuSQUIDS::InteractionsRho(ei,irho) - DRho(ei,irho);
  else
    return nuSQUIDS::InteractionsRho(ei,irho);
}


} // close namespace
