cmake_minimum_required(VERSION 3.2)
enable_testing()
project(afmize)

execute_process(COMMAND git submodule update --init --recursive
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)
include_directories(${PROJECT_SOURCE_DIR}/include)

option(FIND_BOOST "finds boost >= 1.69.0" OFF)

if(BOOST_ROOT OR FIND_BOOST)
    find_package(Boost 1.69.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0/boost/version.hpp")
        message(STATUS "boost 1.72.0 exists.")
        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0)
    else()
        message(STATUS "Boost library not found. Donwloading...")

        file(DOWNLOAD https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz
            ${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0.tar.gz
            EXPECTED_HASH SHA256=c66e88d5786f2ca4dbebb14e06b566fb642a1a6947ad8cc9091f9f445134143f
            STATUS MJOLNIR_DOWNLOAD_BOOST_STATUS)

        # check status
        list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 0 MJOLNIR_DOWNLOAD_BOOST_RESULT)
        if(NOT ${MJOLNIR_DOWNLOAD_BOOST_RESULT} EQUAL 0)
            list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 1 MJOLNIR_DOWNLOAD_BOOST_ERROR)
            message(FATAL_ERROR "failed to download Boost library."
                ${MJOLNIR_DOWNLOAD_BOOST_ERROR})
        endif()

        message(STATUS "downloading completed. Unpacking...")

        execute_process(COMMAND tar xf boost_1_72_0.tar.gz
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0)
        message(STATUS "done.")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
