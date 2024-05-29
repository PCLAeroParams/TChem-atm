#include <TChem_KineticModelData.hpp>
#include <TChem_KineticModelNCAR_ConstData.hpp>
#include <TChem_Util.hpp>

//extern "C" void initialize(const char * filename);
//extern "C" void finalize();

namespace TChem {
struct Driver {
public:
   using real_type_2d_view_host = TChem::real_type_2d_view_host;
   using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
   using device_type      = typename Tines::UseThisDevice<exec_space>::type;

   std::string _chem_file, _therm_file;
   real_type_2d_view_host _state;

   void createStateVector();
   ordinal_type getLengthOfStateVector() const;
   auto getStateVectorSize();
   auto getStateVector();
   void setStateVector(double *array);
   void getStateVectorHost(real_type_2d_const_view_host &view);
 

   // Gases
   TChem::KineticModelData _kmd;
   TChem::KineticModelNCAR_ConstData<host_device_type> _kmcd;   
   void createGasKineticModel(const std::string &chem_file);
   void createGasKineticModelConstData();

   // Aerosols


   // Get sizes
   ordinal_type getNumberOfSpecies();
   std::string getSpeciesName(int *index);


   void doTimestep();
   // Clean up
   void freeAll();
   void freeGasKineticModel();

   // Diagnostics
   void printSpecies();

   // time integration
   real_type_1d_view _t;
   real_type_1d_view _dt;
   real_type_2d_view _tol_time;
   real_type_1d_view _tol_newton;
   real_type_2d_view _fac;

   void setTimeAdvance();

};
}

extern "C" void initialize(const char * filename);
extern "C" void finalize();
extern "C" TChem::ordinal_type TChem_getNumberOfSpecies();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
extern "C" int TChem_getLengthOfStateVector();
extern "C" void TChem_getStateVector(TChem::real_type *array);
extern "C" void TChem_setStateVector(TChem::real_type *array);
extern "C" int TChem_getSpeciesName(int* index, char* result,
     const std::string::size_type buffer_size);
extern "C" void TChem_doTimestep();
extern "C" int TChem_getStateVectorSize();

extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
