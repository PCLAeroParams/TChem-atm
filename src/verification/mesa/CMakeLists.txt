SET(MESA_VALIDATION_DIR ${TCHEM_VERIFICATION_DIR}/mesa)

SET(DRIVER_NAME "mesa_driver")
add_executable(${DRIVER_NAME}.x
               ${DRIVER_NAME}.cpp
               adjust_solid_aerosol.cpp)
TARGET_LINK_LIBRARIES(${DRIVER_NAME}.x verification;${TCHEM_ATM_LINK_LIBRARIES})


# list of test to run
SET(TEST_LIST adjust_solid_aerosol_input_ts_0)
SET(TEST_BASELINE adjust_solid_aerosol_output_ts_0)

SET(DEFAULT_TOL 1e-13)
set(ERROR_THRESHOLDS
    10 # The test is failing because the routine is not invoked. A future PR will invoke it.
   )

foreach(input tol baseline IN ZIP_LISTS TEST_LIST ERROR_THRESHOLDS TEST_BASELINE)
  # add a test to run the skywalker driver
  ADD_TEST(run_${input} ${DRIVER_NAME}.x ${MESA_VALIDATION_DIR}/${input}.yaml)

  # copy the baseline file into place.
   configure_file(
     ${MESA_VALIDATION_DIR}/${baseline}.py
     ${CMAKE_CURRENT_BINARY_DIR}/${input}.py
     COPYONLY
   )
  # add a test to validate TChem-atm's results against the baseline.
  #FIXME: make a version of compare_mam4xx_mam4 for TChem-atm.
  ADD_TEST(validate_${input} python ${TCHEM_VERIFICATION_DIR}/compare_mam4xx_mam4.py
            TChem_${input}.py ${input}.py
            True
            ${DEFAULT_TOL}
            )

endforeach()
