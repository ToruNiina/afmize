cmake_minimum_required(VERSION 3.2)
enable_testing()
project(afmize)

execute_process(COMMAND git submodule update --init --recursive)
include_directories(${PROJECT_SOURCE_DIR})

if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0/boost/version.hpp")
    message(STATUS "boost 1.67.0 exists.")
else()
    message(STATUS "downloading Boost 1.67.0 to ./extlib/ ...")
    execute_process(COMMAND wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
        OUTPUT_QUIET ERROR_QUIET)

    execute_process(COMMAND which sha256sum RESULT_VARIABLE SHA256SUM_EXISTS OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND which shasum    RESULT_VARIABLE SHASUM_EXISTS    OUTPUT_QUIET ERROR_QUIET)
    if(${SHA256SUM_EXISTS} EQUAL "0")
        execute_process(COMMAND sha256sum --check boost_1_67_0_tar_bz2_sha256sum.dat
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
    elseif(${SHASUM_EXISTS} EQUAL "0")
        execute_process(COMMAND shasum --algorithm 256 --check boost_1_67_0_tar_bz2_sha256sum.dat
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
    else()
        message(FATAL_ERROR "neither `sha256sum` nor `shasum` found."
                            "in terms of security, external libraries could not be downloaded."
                            "please install at least one of them.")
    endif()

    if(NOT ${DOWNLOADED_SHA256SUM_RESULT} EQUAL "0")
        message(FATAL_ERROR "boost_1_67_0.tar.bz2 has invalid sha256sum.")
    endif()

    execute_process(COMMAND tar xf boost_1_67_0.tar.bz2
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib" OUTPUT_QUIET ERROR_QUIET)

    message(STATUS "done.")
endif()
include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0)

add_subdirectory(src)
add_subdirectory(test)
