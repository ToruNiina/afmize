add_executable(afmize afmize.cpp)
set_target_properties(afmize PROPERTIES
    COMPILE_FLAGS "-std=c++14 -O2 -DMAVE_NO_SIMD"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

install(TARGETS afmize RUNTIME DESTINATION "bin")
