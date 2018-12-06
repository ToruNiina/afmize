cmake_minimum_required(VERSION 3.2)
enable_testing()
project(afmize)

execute_process(COMMAND git submodule update --init --recursive)
include_directories(${PROJECT_SOURCE_DIR}/extlib/toml)
include_directories(${PROJECT_SOURCE_DIR})

if(BOOST_ROOT)
    find_package(Boost 1.67.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0/boost/version.hpp")
            message(STATUS "boost 1.67.0 exists.")
    else()
        message(STATUS "downloading Boost 1.67.0 to ./extlib/ ...")
        set(BOOST_TAR_BZ2_URL "https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2")

        execute_process(COMMAND which wget RESULT_VARIABLE WGET_EXISTS OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND which curl RESULT_VARIABLE CURL_EXISTS OUTPUT_QUIET ERROR_QUIET)
        if(${WGET_EXISTS} EQUAL "0")
            execute_process(COMMAND wget ${BOOST_TAR_BZ2_URL} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib")
        elseif(${CURL_EXISTS} EQUAL "0")
            execute_process(COMMAND curl ${BOOST_TAR_BZ2_URL} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib")
        else()
            message(FATAL_ERROR "neither wget nor curl found."
                "cannot download boost_X_YY_Z.tar.bz automatically."
                "please install wget or curl, or manually install boost and "
                "provide -DBOOST_ROOT=/path/to/boost when you excute cmake.")
        endif()

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
endif()

add_subdirectory(src)
add_subdirectory(test)
