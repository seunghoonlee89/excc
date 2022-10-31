add_library(cc SHARED ci_to_cc.c ccsd_t.c ecccsd_t.c)

set_target_properties(cc PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})

target_link_libraries(cc ${BLAS_LIBRARIES})