cmake_minimum_required(VERSION 3.2)
enable_testing()
project(afmize)

execute_process(COMMAND git submodule update --init --recursive
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib/toml)
include_directories(${PROJECT_SOURCE_DIR})

if(BOOST_ROOT)
    find_package(Boost 1.69.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0/boost/version.hpp")
        message(STATUS "boost 1.69.0 exists.")
        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0)
    else()
        message(STATUS "Boost library not found. Donwloading...")

        file(DOWNLOAD https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2
            ${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0.tar.bz2
            EXPECTED_HASH SHA256=8f32d4617390d1c2d16f26a27ab60d97807b35440d45891fa340fc2648b04406
            STATUS MJOLNIR_DOWNLOAD_BOOST_STATUS)

        # check status
        list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 0 MJOLNIR_DOWNLOAD_BOOST_RESULT)
        if(NOT ${MJOLNIR_DOWNLOAD_BOOST_RESULT} EQUAL 0)
            list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 1 MJOLNIR_DOWNLOAD_BOOST_ERROR)
            message(FATAL_ERROR "failed to download Boost library."
                ${MJOLNIR_DOWNLOAD_BOOST_ERROR})
        endif()

        message(STATUS "downloading completed. Unpacking...")

        execute_process(COMMAND tar xf boost_1_69_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0)
        message(STATUS "done.")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
