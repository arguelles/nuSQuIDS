#include "nuSQuIDS/nuSQuIDSLV.h"

namespace nusquids {

void nuSQUIDSLV::initialize_LVP_evol_at_each_energy_node() {
    LVP_evol.resize(ne);  // Resize the vector to the number of energy nodes
    for (unsigned int ei = 0; ei < ne; ei++) {
        LVP_evol[ei] = squids::SU_vector(nsun);  // Initialize each SU(N) operator
    }
}

nuSQUIDSLV::nuSQUIDSLV(marray<double, 1> E_range, unsigned int numneu, NeutrinoType NT, bool iinteraction, std::shared_ptr<CrossSectionLibrary> ncs)
    : nuSQUIDS(E_range, numneu, NT, iinteraction, ncs) {
    initialize_LVP_evol_at_each_energy_node();
}

nuSQUIDSLV::nuSQUIDSLV(std::string hdf5_filename, std::string grp, std::shared_ptr<InteractionStructure> int_struct)
    : nuSQUIDS(hdf5_filename, grp, int_struct) {
    // AddToReadHDF5 and AddWriteToHDF5 will be used in constructor of base class
}

nuSQUIDSLV::nuSQUIDSLV(unsigned int numneu, NeutrinoType NT)
    : nuSQUIDS(numneu, NT) {
    initialize_LVP_evol_at_each_energy_node();
}

void nuSQUIDSLV::AddToPreDerive(double x) {
    if (!lv_parameters_set)
        throw std::runtime_error("LV parameters not set");
    for (unsigned int ei = 0; ei < ne; ei++) {
        squids::SU_vector h0 = H0(E_range[ei], 0);
        LVP_evol[ei] = LVP.Evolve(h0, (x - Get_t_initial()));
    }
}

void nuSQUIDSLV::AddToWriteHDF5(hid_t hdf5_loc_id) const {
    double cmutaur = GSL_REAL(c_params.c_mutau);
    double cmutaui = GSL_IMAG(c_params.c_mutau);
    H5LTset_attribute_double(hdf5_loc_id, "c_values", "c_mu_tau_real", &cmutaur, 1);
    H5LTset_attribute_double(hdf5_loc_id, "c_values", "c_mu_tau_imag", &cmutaui, 1);
    double cemur = GSL_REAL(c_params.c_emu);
    double cemui = GSL_IMAG(c_params.c_emu);
    H5LTset_attribute_double(hdf5_loc_id, "c_values", "c_e_mu_real", &cemur, 1);
    H5LTset_attribute_double(hdf5_loc_id, "c_values", "c_e_mu_imag", &cemui, 1);
}

void nuSQUIDSLV::AddToReadHDF5(hid_t hdf5_loc_id) {
    double cmutaur, cmutaui, cemur, cemui;
    H5LTget_attribute_double(hdf5_loc_id, "c_values", "c_mu_tau_real", &cmutaur);
    H5LTget_attribute_double(hdf5_loc_id, "c_values", "c_mu_tau_imag", &cmutaui);
    H5LTget_attribute_double(hdf5_loc_id, "c_values", "c_e_mu_real", &cemur);
    H5LTget_attribute_double(hdf5_loc_id, "c_values", "c_e_mu_imag", &cemui);
    c_params = {gsl_complex_rect(cemur, cemui), gsl_complex_rect(cmutaur, cmutaui)};
    Set_LV_OpMatrix(c_params);
}

squids::SU_vector nuSQUIDSLV::HI(unsigned int ie, unsigned int irho) const {
    squids::SU_vector potential = nuSQUIDS::HI(ie, irho);
    double sign = 1;
    if ((irho == 1 && NT == both) || NT == antineutrino) {
        sign *= -1;
    }
    potential += sign * pow(E_range[ie], n_) * LVP_evol[ie];
    return potential;
}

void nuSQUIDSLV::Set_LV_OpMatrix(LVParameters &lv_params) {
    gsl_matrix_complex *M = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex_set(M, 1, 0, lv_params.c_emu);
    gsl_matrix_complex_set(M, 0, 1, gsl_complex_conjugate(lv_params.c_emu));
    gsl_matrix_complex_set(M, 2, 1, lv_params.c_mutau);
    gsl_matrix_complex_set(M, 1, 2, gsl_complex_conjugate(lv_params.c_mutau));
    LVP = squids::SU_vector(M);
    LVP.RotateToB1(params);
    gsl_matrix_complex_free(M);
    c_params = lv_params;
    lv_parameters_set = true;
}

void nuSQUIDSLV::Set_LV_OpMatrix(gsl_matrix_complex *cmatrix) {
    c_params = {gsl_matrix_complex_get(cmatrix, 1, 0),
                gsl_matrix_complex_get(cmatrix, 2, 1)};
    LVP = squids::SU_vector(cmatrix);
    LVP.RotateToB1(params);
    lv_parameters_set = true;
}

void nuSQUIDSLV::Set_LV_Operator(squids::SU_vector op) {
    LVP = op;
    LVP.RotateToB1(params);
    lv_parameters_set = true;
}

void nuSQUIDSLV::Set_LV_EnergyPower(int n) {
    n_ = n;
}

void nuSQUIDSLV::Set_MixingAngle(unsigned int i, unsigned int j, double angle) {
    nuSQUIDS::Set_MixingAngle(i, j, angle);
    lv_parameters_set = false;
}

void nuSQUIDSLV::Set_CPPhase(unsigned int i, unsigned int j, double angle) {
    nuSQUIDS::Set_CPPhase(i, j, angle);
    lv_parameters_set = false;
}

void nuSQUIDSLV::dump_probabilities() const {
    for (unsigned int ie = 0; ie < ne; ie++) {
        std::cout << E_range[ie] << ' ';
        for (unsigned int flv = 0; flv < numneu; flv++) {
            if (NT == both) {
                std::cout << EvalFlavorAtNode(flv, ie, 0) << ' ';
                std::cout << EvalFlavorAtNode(flv, ie, 1) << ' ';
            } else if (NT == neutrino) {
                std::cout << EvalFlavorAtNode(flv, ie, 0) << ' ';
                std::cout << 0.0 << ' ';
            } else if (NT == antineutrino) {
                std::cout << 0.0 << ' ';
                std::cout << EvalFlavorAtNode(flv, ie, 1) << ' ';
            }
        }
        std::cout << '\n';
    }
}

void nuSQUIDSLV::Set_initial_state(const marray<double, 1> &v, Basis basis) {
    bool lvps = lv_parameters_set;
    nuSQUIDS::Set_initial_state(v, basis);
    lv_parameters_set = lvps;
}

void nuSQUIDSLV::Set_initial_state(const marray<double, 2> &v, Basis basis) {
    bool lvps = lv_parameters_set;
    nuSQUIDS::Set_initial_state(v, basis);
    lv_parameters_set = lvps;
}

void nuSQUIDSLV::Set_initial_state(const marray<double, 3> &v, Basis basis) {
    bool lvps = lv_parameters_set;
    nuSQUIDS::Set_initial_state(v, basis);
    lv_parameters_set = lvps;
}

} // namespace nusquids