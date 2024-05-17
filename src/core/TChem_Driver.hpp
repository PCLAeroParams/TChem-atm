#include <TChem_KineticModelData.hpp>
#include <TChem_KineticModelNCAR_ConstData.hpp>
#include <TChem_Util.hpp>

extern "C" void initialize(const char * filename);
extern "C" void finalize();

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
   auto getStateVector(); //TChem::real_type &view);
   void setStateVector(double *array);

   // 
   TChem::KineticModelData _kmd;
   TChem::KineticModelNCAR_ConstData<host_device_type> _kmcd;   
   void createGasKineticModel(const std::string &chem_file);
   void createGasKineticModelConstData();

   ordinal_type getNumberOfSpecies();
   std::string getSpeciesNames();
   char *getSpeciesName(int *index); //ordinal_type &index);

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

extern "C" TChem::ordinal_type TChem_getNumberOfSpecies();
extern "C" void TChem_getAllStateVectorHost(TChem::real_type *view);
extern "C" const char* TChem_getSpeciesName(TChem::ordinal_type *index);
extern "C" int TChem_getLengthOfStateVector();
extern "C" void TChem_getStateVector(TChem::real_type *array);
extern "C" void TChem_setStateVector(TChem::real_type *array);
extern "C" void TChem_getSpeciesNames(std::string *string);
