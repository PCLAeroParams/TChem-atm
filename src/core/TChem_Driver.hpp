#include <TChem_KineticModelData.hpp>
#include <TChem_KineticModelNCAR_ConstData.hpp>
#include <TChem_AerosolModelData.hpp>
#include <TChem_Util.hpp>

namespace TChem {
struct Driver {
public:
   using real_type_2d_view_host = TChem::real_type_2d_view_host;
   using host_exec_space = Kokkos::DefaultHostExecutionSpace;
   using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
   using device_type = typename Tines::UseThisDevice<exec_space>::type;

   // Input files
   std::string _chem_file, _aero_file;

   ordinal_type _team_size;
   ordinal_type _vector_size;

   // State vector
   ordinal_type _nBatch;
   void setBatchSize(ordinal_type nBatch);
   real_type_2d_view_host _state;
   void createStateVector(ordinal_type nBatch);
   ordinal_type getLengthOfStateVector() const;
   auto getStateVectorSize();
   auto getStateVector(const ordinal_type iBatch);
   void setStateVector(double *array, const ordinal_type iBatch);
   void getStateVectorHost(real_type_2d_const_view_host &view);

   // Gas variables
   TChem::KineticModelData _kmd;
   TChem::KineticModelNCAR_ConstData<host_device_type> _kmcd_host;
   TChem::KineticModelNCAR_ConstData<device_type> _kmcd_device;
   void createGasKineticModel(const std::string &chem_file);
   void createGasKineticModelConstData();

   // Return number of gas species
   ordinal_type getNumberOfSpecies();
   // Return gas species name
   std::string getSpeciesName(int *index);

   // Aerosol variables
   TChem::AerosolModelData _amd;
   TChem::AerosolModel_ConstData<host_device_type> _amcd_host;
   TChem::AerosolModel_ConstData<device_type> _amcd_device;
   void createAerosolModel(const std::string &aero_file);
   void createAerosolModelConstData();

   // Return number of aerosol species
   ordinal_type getNumberOfAerosolSpecies();

   // Return aerosol species information
   real_type getAerosolSpeciesDensity(int *index);
   real_type getAerosolSpeciesMW(int *index);
   real_type getAerosolSpeciesKappa(int *index);

   // Integrate a single time step
   void doTimestep(const double del_t);

   // Time integration information
   real_type _atol_newton; // Absolute tolerance used in Newton solver
   real_type _rtol_newton; // Relative tolerance used in Newton solver
   real_type _dtmin; // Minimum time step size (s).
   real_type _atol_time; // Absolute tolerance used for adaptive time stepping
   real_type _rtol_time; // Relative tolerance used for adaptive time stepping
   ordinal_type _max_num_newton_iterations; // Maximum number of Newton iterations
   ordinal_type _max_num_time_iterations; // Maximum number of time iterations
   ordinal_type _num_time_iterations_per_interval; // 
   ordinal_type _jacobian_interval; // Interval for evaluating Jacobians

   // Read in time integration information
   void createNumerics(const std::string &numerics_file);

   // Clean up
   void freeAll();
   void freeGasKineticModel();
   void freeAerosolModel();

   // Diagnostics

};
}

extern "C" void initialize(const char* gasFile, const char* aeroFile,
                           const char* numericsFile, const TChem::ordinal_type nBatch);
extern "C" void finalize();
extern "C" TChem::ordinal_type TChem_getNumberOfSpecies();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
extern "C" int TChem_getLengthOfStateVector();
extern "C" void TChem_getStateVector(TChem::real_type *array, const TChem::ordinal_type iBatch);
extern "C" void TChem_setStateVector(TChem::real_type *array, const TChem::ordinal_type iBatch);
extern "C" int TChem_getSpeciesName(int* index, char* result,
                                    const std::string::size_type buffer_size);
extern "C" void TChem_doTimestep(const double &del_t);
extern "C" int TChem_getStateVectorSize();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);

extern "C" TChem::ordinal_type TChem_getNumberOfAerosolSpecies();
extern "C" int TChem_getAerosolSpeciesName(int* index, char* result,
                                    const std::string::size_type buffer_size);
extern "C" double TChem_getAerosolSpeciesDensity(int* index);
extern "C" double TChem_getAerosolSpeciesMW(int* index);
extern "C" double TChem_getAerosolSpeciesKappa(int* index);
