cmake_minimum_required(VERSION 3.2)
enable_testing()
project(afmize)

execute_process(COMMAND git submodule update --init --recursive
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)
include_directories(${PROJECT_SOURCE_DIR}/include)


find_package(Boost 1.69.0)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0/boost/version.hpp")
        message(STATUS "boost 1.72.0 exists.")
        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_72_0)
    else()
        message(STATUS "Boost library not found. Donwloading...")

        file(DOWNLOAD https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz
            ${PROJECT_SOURCE_DIR}/extlib/boost_1_77_0.tar.gz
            EXPECTED_HASH SHA256=5347464af5b14ac54bb945dc68f1dd7c56f0dad7262816b956138fc53bcc0131
            STATUS AFMIZE_DOWNLOAD_BOOST_STATUS)

        # check status
        list(GET AFMIZE_DOWNLOAD_BOOST_STATUS 0 AFMIZE_DOWNLOAD_BOOST_RESULT)
        if(NOT ${AFMIZE_DOWNLOAD_BOOST_RESULT} EQUAL 0)
            list(GET AFMIZE_DOWNLOAD_BOOST_STATUS 1 AFMIZE_DOWNLOAD_BOOST_ERROR)
            message(FATAL_ERROR "failed to download Boost library."
                ${AFMIZE_DOWNLOAD_BOOST_ERROR})
        endif()

        message(STATUS "downloading completed. Unpacking...")

        execute_process(COMMAND tar xf boost_1_77_0.tar.gz
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_77_0)
        message(STATUS "done.")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
