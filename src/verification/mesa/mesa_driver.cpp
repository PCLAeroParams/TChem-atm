#include <iostream>
#include <verification.hpp>

#if defined(TCHEM_ATM_ENABLE_SKYWALKER)
 #include "skywalker.hpp"
 using namespace skywalker;
 using namespace tchem;
#endif
void usage() {
  std::cerr << "mesa_driver: a Skywalker driver for validating the "
               "mesa routines."
            << std::endl;
  std::cerr << "mesa_driver: usage:" << std::endl;
  std::cerr << "mesa_driver <input.yaml>" << std::endl;
  exit(0);
}
// void adjust_solid_aerosol(Ensemble *ensemble);

int main(int argc, char **argv) {
  if (argc == 1) {
    usage();
  }

  verification::initialize(argc, argv);
  std::string input_file = argv[1];
  std::string output_file = verification::output_name(input_file);
  std::cout << argv[0] << ": reading " << input_file << std::endl;

  // Load the ensemble. Any error encountered is fatal.
  Ensemble *ensemble = skywalker::load_ensemble(input_file, "mam4xx");

  // the settings.
  Settings settings = ensemble->settings();
  if (!settings.has("function")) {
    std::cerr << "No function specified in mam4xx.settings!" << std::endl;
    exit(1);
  }

  // Dispatch to the requested function.
  auto func_name = settings.get("function");
  try {
#if 0
    if (func_name == "adjust_solid_aerosol") {
      adjust_solid_aerosol(ensemble);
    } else {
      std::cerr << "Error: Function name '" << func_name
                << "' does not have an implemented test!" << std::endl;
      exit(1);
    }
#endif
  } catch (std::exception &e) {
    std::cerr << argv[0] << ": Error: " << e.what() << std::endl;
  }

  // Write out a Python module.
  std::cout << argv[0] << ": writing " << output_file << std::endl;
  ensemble->write(output_file);

  // Clean up.
  delete ensemble;
  verification::finalize();
}
