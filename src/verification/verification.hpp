#ifndef TCHEM_ATM_VERIFICATION_HPP
#define TCHEM_ATM_VERIFICATION_HPP

#include "TChem.hpp"
namespace tchem {
namespace verification {

/// Call this function to initialize a validation driver.
void initialize(int argc, char **argv);

/// Call this function to finalize a validation driver.
void finalize();

std::string output_name(const std::string &input_file);
}// namespace verification
} // namespace tchem

#endif
