
using namespace nusquids;


class nuSQUIDSNSI: public nuSQUIDS {
  private:
    squids::SU_vector NSI;
    std::vector<squids::SU_vector> NSI_evol;
    std::unique_ptr<double[]> hiBuffer;
    double HI_prefactor;
    // nsi parameters
    double epsilon_mutau;

    void AddToPreDerive(double x){
      for(int ei = 0; ei < ne; ei++){
        // asumming same hamiltonian for neutrinos/antineutrinos
        //SU_vector h0 = H0(E_range[ei],0);
        //NSI_evol[ei] = NSI.Evolve(h0,(x-Get_t_initial()));
        NSI_evol[ei] = NSI.Evolve(H0_array[ei],(x-Get_t_initial()));
      }
    }

    void AddToReadHDF5(hid_t hdf5_loc_id){
      // here we read the new parameters now saved in the HDF5 file
      hid_t nsi = H5Gopen(hdf5_loc_id, "nsi", H5P_DEFAULT);
      H5LTget_attribute_double(hdf5_loc_id,"nsi","mu_tau" ,&epsilon_mutau);
      H5Gclose(nsi);
    }

    void AddToWriteHDF5(hid_t hdf5_loc_id) const {
      // here we write the new parameters to be saved in the HDF5 file
      H5Gcreate(hdf5_loc_id, "nsi", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5LTset_attribute_double(hdf5_loc_id, "nsi","mu_tau",&epsilon_mutau, 1);
    }

    squids::SU_vector HI(unsigned int ei,unsigned int index_rho) const{
      double CC = HI_prefactor*body->density(*track)*body->ye(*track);

      // // construct potential in flavor basis
      squids::SU_vector potential(nsun,hiBuffer.get());

      potential = (3.0*CC)*NSI_evol[ei];

      if ((index_rho == 0 and NT==both) or NT==neutrino){
          // neutrino potential
          return nuSQUIDS::HI(ei,index_rho) + potential;
      } else if ((index_rho == 1 and NT==both) or NT==antineutrino){
          // antineutrino potential
          return nuSQUIDS::HI(ei,index_rho) + (-1.0)*std::move(potential);
      } else{
          throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
      }
    }
  public:
  //void constructor, is needed for some wierd reason beyon my capability 
  //of understanding c++ but probably Chris knows
//  nuSQUIDSNSI(){}
  //Constructor
  nuSQUIDSNSI(double epsilon_mutau, marray<double,1> Erange,unsigned int numneu, NeutrinoType NT,
              bool iinteraction,double th01=0.563942, double th02=0.154085,double th12=0.785398):
                      nuSQUIDS(Erange,numneu,NT,iinteraction),
				      hiBuffer(new double[nsun*nsun]),epsilon_mutau(epsilon_mutau)
  {
    assert(numneu == 3);
    // defining a complex matrix M which will contain our flavor
    // violating flavor structure.
    gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
    gsl_complex c {{ epsilon_mutau , 0.0 }};
    gsl_matrix_complex_set(M,2,1,c);
    gsl_matrix_complex_set(M,1,2,gsl_complex_conjugate(c));
    
    NSI = squids::SU_vector(M);

    //setting the mising angle for the flavro-mass basis, 
    //they are measured by the neutrino oscillation experiments.    
    Set_MixingAngle(0,1,th01);
    Set_MixingAngle(0,2,th02);
    Set_MixingAngle(1,2,th12);
    
    // rotate to mass reprentation
    NSI.RotateToB1(params);
    NSI_evol.resize(ne);
    for(int ei = 0; ei < ne; ei++){
      NSI_evol[ei] = squids::SU_vector(nsun);
    }
    gsl_matrix_complex_free(M);
    
    HI_prefactor = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);
  }
  
  void Set_mutau(double eps){
    gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
    gsl_complex c {{ epsilon_mutau , 0.0 }};
    gsl_matrix_complex_set(M,2,1,c);
    gsl_matrix_complex_set(M,1,2,gsl_complex_conjugate(c));
    NSI = squids::SU_vector(M);    
    NSI.RotateToB1(params);
    gsl_matrix_complex_free(M);
  }
  
  
};
