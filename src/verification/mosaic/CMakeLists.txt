SET(MOSAIC_VALIDATION_DIR ${TCHEM_VERIFICATION_DIR}/mosaic)

SET(DRIVER_NAME "mosaic_driver")
add_executable(${DRIVER_NAME}.x
               ${DRIVER_NAME}.cpp
               adjust_liquid_aerosol.cpp
               adjust_solid_aerosol.cpp
               do_full_deliquescence.cpp
               calculate_XT.cpp
               fnlog_gamZ.cpp
               fn_Keq.cpp
               fn_Po.cpp
               molality_0.cpp
               bin_molality.cpp
               bin_molality_60.cpp
               MTEM_compute_log_gamZ.cpp
 )



TARGET_LINK_LIBRARIES(${DRIVER_NAME}.x verification;${TCHEM_ATM_LINK_LIBRARIES})


# list of test to run
SET(TEST_LIST
     adjust_liquid_aerosol_input_ts_0
     adjust_solid_aerosol_input_ts_0
     do_full_deliquescence_input_ts_0
     calculate_XT_input_ts_0
     fnlog_gamZ_input_ts_0
     fn_Keq_input_ts_0
     fn_Po_input_ts_0
     molality_0_input_ts_0
     bin_molality_input_ts_0
     bin_molality_60_input_ts_0
     MTEM_compute_log_gamZ_input_ts_0
     aerosol_water_up_input_ts_0
)

SET(TEST_BASELINE
    adjust_liquid_aerosol_output_ts_0
    adjust_solid_aerosol_output_ts_0
    do_full_deliquescence_output_ts_0
    calculate_XT_output_ts_0
    fnlog_gamZ_output_ts_0
    fn_Keq_output_ts_0
    fn_Po_output_ts_0
    molality_0_output_ts_0
    bin_molality_output_ts_0
    bin_molality_60_output_ts_0
    MTEM_compute_log_gamZ_output_ts_0
    aerosol_water_up_output_ts_0
)

SET(DEFAULT_TOL 1e-9)
set(ERROR_THRESHOLDS
    ${DEFAULT_TOL} # adjust_liquid_aerosol_output_ts_0
    ${DEFAULT_TOL} # adjust_solid_aerosol_output_ts_0
    ${DEFAULT_TOL} # do_full_deliquescence_output_ts_0
    ${DEFAULT_TOL} # calculate_XT_output_ts_0
    1e-6           # fnlog_gamZ_output_ts_0
    1e-6           # fn_Keq_output_ts_0
    1e-8           # fn_Po_output_ts_0
    1e-5           # molality_0_output_ts_0
    1e-5           # bin_molality_output_ts_0
    1e-6           # bin_molality_60_output_ts_0
    1.8e-4         # MTEM_compute_log_gamZ
    ${DEFAULT_TOL} # aerosol_water_up_output_ts_0
)

foreach(input tol baseline IN ZIP_LISTS TEST_LIST ERROR_THRESHOLDS TEST_BASELINE)
  # add a test to run the skywalker driver
  ADD_TEST(run_${input} ${DRIVER_NAME}.x ${MOSAIC_VALIDATION_DIR}/${input}.yaml)

  # copy the baseline file into place.
   configure_file(
     ${MOSAIC_VALIDATION_DIR}/${baseline}.py
     ${CMAKE_CURRENT_BINARY_DIR}/${input}.py
     COPYONLY
   )
  # add a test to validate TChem-atm's results against the baseline.
  ADD_TEST(validate_${input} python3 ${TCHEM_VERIFICATION_DIR}/file_difference_checker.py
            TChem_${input}.py ${input}.py
            True 
            ${tol}
            )

endforeach()
