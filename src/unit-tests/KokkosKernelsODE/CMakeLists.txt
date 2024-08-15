# Copy the 'references' directory to the current binary directory
FILE(COPY references
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Set the test name and executable name
SET(TEST_NAME AerosolChemistry_KokkosKernels)
SET(EXE_NAME ../../examples/TChem_${TEST_NAME}.x)

# Set the input path for the example run
SET(INPUT_PATH ../../examples/runs/atmospheric_chemistry/Simpol)

# Set the validation script path
SET(VALIDATION_SCRIPT ../scripts/compare_tchem.py)

# Threshold error for the passing test
SET(PASSING_ERROR 1e-8)

# Set the output file name
SET(OUTPUT_NAME full_gas.dat)

# Add a test to execute the example code
ADD_TEST(run_${TEST_NAME} ${EXE_NAME} 
    --chemfile=${INPUT_PATH}/config_gas.yaml 
    --aerofile=${INPUT_PATH}/test_SIMPOL_phase_transfer.yaml 
    --inputfile_particles=${INPUT_PATH}/scenario_conditions_particle.yaml 
    --outputfile=${OUTPUT_NAME} 
    --dtmin=1 
    --tend=5)

# Add a test to validate the output data with the reference data
ADD_TEST(validate_${TEST_NAME} python ${VALIDATION_SCRIPT} 
    -r_file=${OUTPUT_NAME}
    -t_file=references/${OUTPUT_NAME}
    -check_norms=True
    =error=${PASSING_ERROR})
