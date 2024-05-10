#include <TChem_KineticModelData.hpp>
#include <TChem_KineticModelNCAR_ConstData.hpp>

extern "C" void normal_vec(int n, double x[]);
extern "C" void call_something(int n);
extern "C" void do_something(int n);
extern "C" void initialize_kokkos(const char * filename);
extern "C" void finalize_kokkos();

namespace TChem {
struct Driver {
public:
   using real_type_2d_view_host = TChem::real_type_2d_view_host;
   using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
   using device_type      = typename Tines::UseThisDevice<exec_space>::type;

   std::string _chem_file, _therm_file;
   TChem::KineticModelData _kmd;
   TChem::KineticModelNCAR_ConstData<host_device_type> _kmcd;
   real_type_2d_view_host state_host;

   void createGasKineticModel(const std::string &chem_file);
   void createGasKineticModelConstData();

   void freeAll();
   void freeGasKineticModel();

   void printSpecies();

};
}

extern "C" int TChem_getNumberOfSpecies();
