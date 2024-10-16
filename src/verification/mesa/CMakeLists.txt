SET(DRIVER_NAME "mesa_driver")
add_executable(${DRIVER_NAME}.x
               ${DRIVER_NAME}.cpp
               adjust_solid_aerosol.cpp)
TARGET_LINK_LIBRARIES(${DRIVER_NAME}.x verification;${TCHEM_ATM_LINK_LIBRARIES})
