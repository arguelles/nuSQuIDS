
using namespace nusquids;


class nuSQUIDSNSI: public nuSQUIDS {
  private:
    squids::SU_vector NSI;
    squids::SU_vector SI;
    squids::SU_vector custom_HI;
    std::vector<squids::SU_vector> NSI_evol;
    std::vector<squids::SU_vector> SI_evol;
    std::unique_ptr<double[]> hiBuffer;
    double HI_prefactor;
    double DMP;
    void AddToPreDerive(double x){
      for(int ei = 0; ei < ne; ei++){
        //asumming same hamiltonian for neutrinos/antineutrinos
        //SU_vector h0 = H0(E_range[ei],0);
        //NSI_evol[ei] = NSI.Evolve(h0,(x-Get_t_initial()));
        NSI_evol[ei] = NSI.Evolve(H0_array[ei],(x-Get_t_initial()));
        SI_evol[ei] = SI.Evolve(H0_array[ei],(x-Get_t_initial()));
      }
    }
    
    
    squids::SU_vector HI(unsigned int ei,unsigned int index_rho) const{
      squids::SU_vector potentialNSI(nsun,hiBuffer.get());
      squids::SU_vector potentialSI(nsun,hiBuffer.get());

      double CC = HI_prefactor*current_density*current_ye;

      potentialSI = (-1.0)*CC*SI_evol[ei];
      potentialNSI = NSI_evol[ei];



      if ((index_rho == 0 and NT==both) or NT==neutrino){
          // neutrino potential
          return 0*potentialSI + potentialNSI;
      } else if ((index_rho == 1 and NT==both) or NT==antineutrino){
          // antineutrino potential
          return 0*potentialSI + (-1.0)*std::move(potentialNSI);
      } else{
          throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
      }
    }
  public:
  nuSQUIDSNSI(marray<double,1> Erange,int numneu, NeutrinoType NT,
	      bool iinteraction,double th04=0.00305128, double th15=0.003,
	      double th26=0.00295122,double th34=0.3, double th45=0.563942,
	      double th46=0.154085,double th56=0.785398) : nuSQUIDS(Erange,numneu,NT,iinteraction),
				      hiBuffer(new double[nsun*nsun])
  {    
    DMP = -4.16667e-7*2.65/2.65;
    // defining a complex matrix M which will contain our BSM Potential
    gsl_matrix_complex * M = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_complex c {{ DMP , 0.0 }};
    gsl_matrix_complex_set(M,0,0,c);
    gsl_matrix_complex_set(M,1,1,c);
    gsl_matrix_complex_set(M,2,2,c); 
    
    NSI = squids::SU_vector(M);
 
    HI_prefactor = params.sqrt2*params.GF*params.Na*pow(params.cm,-3); 
    gsl_matrix_complex * MSM = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_complex d {{ 1 , 0.0 }};
    gsl_matrix_complex_set(MSM,5,5,d);
    gsl_matrix_complex_set(MSM,6,6,d);
    SI = squids::SU_vector(MSM);
 
 
    Set_MixingAngle(0,4,th04);
    Set_MixingAngle(1,5,th15);
    Set_MixingAngle(2,6,th26);
    Set_MixingAngle(3,4,th34);
    Set_MixingAngle(4,5,th45);
    Set_MixingAngle(4,6,th46);
    Set_MixingAngle(5,6,th56);
    
    // rotate to mass reprentation
    NSI.RotateToB1(params);
    SI.RotateToB1(params);
    NSI_evol.resize(ne);
    SI_evol.resize(ne);
    for(int ei = 0; ei < ne; ei++){
      NSI_evol[ei] = squids::SU_vector(nsun);
      SI_evol[ei] = squids::SU_vector(nsun);
    }
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(MSM);

  }

  
  
};
