#include <TChem_KineticModelData.hpp>
#include <TChem_KineticModelNCAR_ConstData.hpp>
#include <TChem_Util.hpp>

namespace TChem {
struct Driver {
public:
   using real_type_2d_view_host = TChem::real_type_2d_view_host;
   using host_exec_space = Kokkos::DefaultHostExecutionSpace;
   using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
   using device_type = typename Tines::UseThisDevice<exec_space>::type;

   // Input files
   std::string _chem_file, _therm_file;

   // State vector
   real_type_2d_view_host _state;
   void createStateVector();
   ordinal_type getLengthOfStateVector() const;
   auto getStateVectorSize();
   auto getStateVector();
   void setStateVector(double *array);
   void getStateVectorHost(real_type_2d_const_view_host &view);

   // Gas variables
   TChem::KineticModelData _kmd;
   TChem::KineticModelNCAR_ConstData<host_device_type> _kmcd_host;
   TChem::KineticModelNCAR_ConstData<device_type> _kmcd_device;
   void createGasKineticModel(const std::string &chem_file);
   void createGasKineticModelConstData();

   // Aerosols

   // Return number of gas species
   ordinal_type getNumberOfSpecies();
   // Return gasspecies name
   std::string getSpeciesName(int *index);

   // Integrate a single time step
   void doTimestep(const double del_t);

   // Time integration information
   real_type _atol_newton;
   real_type _rtol_newton;
   real_type _dtmin;
   real_type _atol_time;
   real_type _rtol_time;
   ordinal_type _max_num_newton_iterations;
   ordinal_type _num_time_iterations_per_interval;
   ordinal_type _jacobian_interval;

   // Read in time integration information
   void createNumerics(const std::string &numerics_file);

   // Clean up
   void freeAll();
   void freeGasKineticModel();

   // Diagnostics

};
}

extern "C" void initialize(const char* gasFile, const char* aeroFile,
                           const char* numericsFile);
extern "C" void finalize();
extern "C" TChem::ordinal_type TChem_getNumberOfSpecies();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
extern "C" int TChem_getLengthOfStateVector();
extern "C" void TChem_getStateVector(TChem::real_type *array);
extern "C" void TChem_setStateVector(TChem::real_type *array);
extern "C" int TChem_getSpeciesName(int* index, char* result,
                                    const std::string::size_type buffer_size);
extern "C" void TChem_doTimestep(const double &del_t);
extern "C" int TChem_getStateVectorSize();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
