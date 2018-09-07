add_executable(afmize afmize.cpp)
# target_link_libraries(afmize )
set_target_properties(afmize PROPERTIES COMPILE_FLAGS "-std=c++14 -O2")
