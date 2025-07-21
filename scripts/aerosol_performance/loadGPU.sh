
# Uncomment to use the Kokkos tools NVTX connector when profiling
#export KOKKOS_TOOLS_LIBS=${KOKKOS_TOOLS_INSTALL_DIR}/lib64/libkp_nvtx_connector.so

# MemoryEvents - tracks timeline of alloc and dealloc events and mem usage
export KOKKOS_TOOLS_LIBS=${KOKKOS_TOOLS_INSTALL_DIR}/lib64/libkp_memory_events.so

TCHEM_INSTALL_PATH=${TCHEM_DEVICE_INSTALL_PATH}
exec=${TCHEM_INSTALL_PATH}/examples/${experiment_exe_root}.x
machine_name=nid

