SET(MESA_VALIDATION_DIR ${TCHEM_VERIFICATION_DIR}/mesa)

SET(DRIVER_NAME "mesa_driver")
add_executable(${DRIVER_NAME}.x
               ${DRIVER_NAME}.cpp
               adjust_solid_aerosol.cpp)
TARGET_LINK_LIBRARIES(${DRIVER_NAME}.x verification;${TCHEM_ATM_LINK_LIBRARIES})


# list of test to run
SET(TEST_LIST adjust_solid_aerosol_input_ts_0)

SET(DEFAULT_TOL 1e-13)
set(ERROR_THRESHOLDS
    ${DEFAULT_TOL}
   )

foreach(input tol IN ZIP_LISTS TEST_LIST ERROR_THRESHOLDS)
  # copy the baseline file into place.

  # add a test to run the skywalker driver
  ADD_TEST(run_${input} ${DRIVER_NAME}.x ${MESA_VALIDATION_DIR}/${input}.yaml)
endforeach()
