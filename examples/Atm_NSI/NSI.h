
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
      double ye = body->ye(*track);
      double density = body->density(*track);

      double CC = HI_prefactor*density*ye;
      double NC;

      if (ye < 1.0e-10){
        NC = HI_prefactor*density;
      }
      else {
        NC = CC*(-0.5*(1.0-ye)/ye);
      }

      // construct potential in flavor basis
      squids::SU_vector potential(nsun,hiBuffer.get());
      potential = (CC+NC)*evol_b1_proj[index_rho][0][ei];
      potential += (NC)*(evol_b1_proj[index_rho][1][ei]);
      potential += (NC)*(evol_b1_proj[index_rho][2][ei]);
      // plus sign so that the NSI potential has the same sign as the VCC potential
      // and the factor of 3 comes from average n_n/n_e at Earth.
      potential += (3.0*CC)*NSI_evol[ei];

      if ((index_rho == 0 and NT==both) or NT==neutrino){
          // neutrino potential
          return potential;
      } else if ((index_rho == 1 and NT==both) or NT==antineutrino){
          // antineutrino potential
          return (-1.0)*std::move(potential);
      } else{
          throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
      }
    }
public:
  nuSQUIDSNSI(){}
  nuSQUIDSNSI(double epsilon_mutau,double Emin,double Emax,int Esize,int numneu, NeutrinoType NT,
	      bool elogscale,bool iinteraction) : nuSQUIDS(Emin,Emax,Esize,numneu,NT,elogscale,iinteraction),
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
    
    Set_MixingAngle(0,1,0.563942);
    Set_MixingAngle(0,2,0.154085);
    Set_MixingAngle(1,2,0.785398);
    
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
