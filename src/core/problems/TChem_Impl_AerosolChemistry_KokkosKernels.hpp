/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_IMPL_AEROSOL_CHEMISTRY_KOKKOSKERNELS_HPP__
#define __TCHEM_IMPL_AEROSOL_CHEMISTRY_KOKKOSKERNELS_HPP__

#include "Tines_Internal.hpp"

#include "TChem_KineticModelData.hpp"
#include "TChem_AerosolModelData.hpp"
#include "TChem_Impl_AerosolChemistry_Problem.hpp"
#include "KokkosODE_BDF.hpp"

namespace TChem {
namespace Impl {

template<typename ValueType, typename DeviceType>
struct StiffChemistry {  

  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<value_type, device_type>;

  // FIXME: let's try first Serial and host_device_type
  using host_device_type = Tines::UseThisDevice<host_exec_space>::type;
  // using problem_type = AerosolChemistry_Problem<realtype,host_device_type>;
  problem_type _problem;
  int neqs;


  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/,
                                         const double /*dt*/,
                                         const vec_type1& y,
                                         const vec_type2& f) const {
   
  // FIXME: use memger in device. 
  const auto member = Tines::HostSerialTeamMember();  
    _problem.computeFunction(member,
                             y,
                             f);
  }

    template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/,
                                         const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
  
  // FIXME: use memger in device. 
  const auto member = Tines::HostSerialTeamMember();  
  _problem.computeJacobian(member,
             y,
             jac);
              
  }


};

} // Impl
}// TChem

#endif    