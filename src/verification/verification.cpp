#include "verification.hpp"

namespace tchem {
namespace verification {
void initialize(int argc, char **argv) { Kokkos::initialize(argc, argv); }

void finalize() {
  finalize();
  Kokkos::finalize();
  }

std::string output_name(const std::string &input_file) {
  std::string output_file;
  size_t slash = input_file.find_last_of('/');
  size_t dot = input_file.find_last_of('.');
  if ((dot == std::string::npos) and (slash == std::string::npos)) {
    dot = input_file.length();
  }
  if (slash == std::string::npos) {
    slash = 0;
  } else {
    slash += 1;
    dot -= slash;
  }
  return std::string("TChem_") + input_file.substr(slash, dot) +
         std::string(".py");
}

void convert_1d_vector_to_2d_view_device(const std::vector<real_type> &var_std,
                                         const real_type_2d_view &var_device) {
  auto host = Kokkos::create_mirror_view(var_device);
  int count = 0;
  for (int d2 = 0; d2 < var_device.extent(1); ++d2) {
    for (int d1 = 0; d1 < var_device.extent(0); ++d1) {
      host(d1, d2) = var_std[count];
      count++;
    }
  }
  Kokkos::deep_copy(var_device, host);
}

void convert_2d_view_device_to_1d_vector(const real_type_2d_view &var_device,
                                         std::vector<real_type> &var_std) {
  auto host = Kokkos::create_mirror_view(var_device);
  Kokkos::deep_copy(host, var_device);
  int count = 0;
  for (int d2 = 0; d2 < var_device.extent(1); ++d2) {
    for (int d1 = 0; d1 < var_device.extent(0); ++d1) {
      var_std[count] = host(d1, d2);
      count++;
    }
  }
}

}   // namespace verification
} // namespace tchem
