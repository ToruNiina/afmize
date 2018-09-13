add_executable(afmize afmize.cpp)
set_target_properties(afmize PROPERTIES
    COMPILE_FLAGS "-std=c++14 -O2"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
