
using namespace nusquids;


class nuSQUIDSNSI: public nuSQUIDS {
  private:
    squids::SU_vector NSI;
    std::vector<squids::SU_vector> NSI_evol;
    std::unique_ptr<double[]> hiBuffer;
    // nsi parameters
    double DMP;

    void AddToPreDerive(double x){
      for(int ei = 0; ei < ne; ei++){
        // asumming same hamiltonian for neutrinos/antineutrinos
        //SU_vector h0 = H0(E_range[ei],0);
        //NSI_evol[ei] = NSI.Evolve(h0,(x-Get_t_initial()));
        NSI_evol[ei] = NSI.Evolve(H0_array[ei],(x-Get_t_initial()));
      }
    }

    squids::SU_vector HI(unsigned int ei,unsigned int index_rho) const{
      squids::SU_vector potential(nsun,hiBuffer.get());

      potential = NSI_evol[ei];

      if ((index_rho == 0 and NT==both) or NT==neutrino){
          // neutrino potential
          return 0*nuSQUIDS::HI(ei,index_rho) + potential;
      } else if ((index_rho == 1 and NT==both) or NT==antineutrino){
          // antineutrino potential
          return 0*nuSQUIDS::HI(ei,index_rho) + (-1.0)*std::move(potential);
      } else{
          throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
      }
    }
  public:
  nuSQUIDSNSI(double DMP, marray<double,1> Erange,int numneu, NeutrinoType NT,
	      bool iinteraction,double th04=0.00305128, double th15=0.003,
	      double th26=0.00295122,double th34=0.3, double th45=0.563942,
	      double th46=0.154085,double th56=0.785398) : nuSQUIDS(Erange,numneu,NT,iinteraction),
				      hiBuffer(new double[nsun*nsun])
  {
    // defining a complex matrix M which will contain our flavor
    // violating flavor structure.
    gsl_matrix_complex * M = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_complex c {{ DMP , 0.0 }};
    gsl_matrix_complex_set(M,0,0,c);
    gsl_matrix_complex_set(M,1,1,c);
    gsl_matrix_complex_set(M,2,2,c);
    
    NSI = squids::SU_vector(M);
    
    Set_MixingAngle(0,4,th04);
    Set_MixingAngle(1,5,th15);
    Set_MixingAngle(2,6,th26);
    Set_MixingAngle(3,4,th34);
    Set_MixingAngle(4,5,th45);
    Set_MixingAngle(4,6,th46);
    Set_MixingAngle(5,6,th56);
    
    // rotate to mass reprentation
    NSI.RotateToB1(params);
    NSI_evol.resize(ne);
    for(int ei = 0; ei < ne; ei++){
      NSI_evol[ei] = squids::SU_vector(nsun);
    }
    gsl_matrix_complex_free(M);
    
  }

  
  
};
