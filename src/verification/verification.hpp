#ifndef TCHEM_ATM_VERIFICATION_HPP
#define TCHEM_ATM_VERIFICATION_HPP

#include "TChem.hpp"

namespace tchem {
using real_type = TChem::real_type;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_1d_view_host = TChem::real_type_1d_view_host;
namespace verification {

/// Call this function to initialize a validation driver.
void initialize(int argc, char **argv);

/// Call this function to finalize a validation driver.
void finalize();

std::string output_name(const std::string &input_file);

void convert_1d_vector_to_2d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_2d_view &var_device);

void convert_2d_view_device_to_1d_vector(const real_type_2d_view &var_device,
                                         std::vector<real_type> &var_std);

void convert_1d_vector_to_1d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_1d_view &var_device);
}// namespace verification
} // namespace tchem

#endif
