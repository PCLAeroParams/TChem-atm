#ifndef __TCHEM_REACTION_TYPES_HPP__
#define __TCHEM_REACTION_TYPES_HPP__

#include "TChem_Util.hpp"
namespace TChem {


struct SIMPOL_PhaseTransferType{
  // index in state vector
  ordinal_type aerosol_sp_index;
  // mass transfer parameters
  real_type B1;
  real_type B2;
  real_type B3;
  real_type B4;
  // gas species info
  ordinal_type gas_sp_index;
  real_type diffusion_coeff;
  real_type N_star;
  bool compute_alpha;
  real_type molecular_weight;
};

using simpol_phase_transfer_type_1d_dual_view =
 Tines::value_type_1d_dual_view<SIMPOL_PhaseTransferType, exec_space>;

struct ToyProcessType {
  // kinetic parameters
  ordinal_type A;
  ordinal_type B;

  // aerosol index
  ordinal_type aerosol_sp_index;

  // gas species info
  ordinal_type gas_sp_index;
};

using toy_process_1d_dual_view =
 Tines::value_type_1d_dual_view<ToyProcessType, exec_space>;

} // namespace TChem

#endif
