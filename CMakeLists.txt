# MIT License
#
# Copyright (c) 2021 Mechatronics and Haptic Interfaces Lab - Rice University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# Author(s): Nathan Dunkelberger (nbd2@rice.edu)

cmake_minimum_required(VERSION 3.13)
project(mahi_mpc VERSION 0.1.0 LANGUAGES CXX)

# Options

option(WITH_MAHI_MPC "Turn ON to build example executable(s)" ON)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
option(MAHI_MPC_EXAMPLES "Turn ON to build example executable(s)" ON)
else()
option(MAHI_MPC_EXAMPLES "Turn ON to build example executable(s)" OFF)
endif()

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS "-openmp")

set(IPOPT_BASE_FILEPATH "C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md")
set(IPOPT_LIBRARIES "${IPOPT_BASE_FILEPATH}/lib/ipopt.dll.lib")
set(IPOPT_INCLUDE_DIRS "${IPOPT_BASE_FILEPATH}/include/coin-or")
set(Ipopt_DIR "${IPOPT_BASE_FILEPATH}/bin")
set(WITH_IPOPT ON CACHE BOOL "IPOPT will be used" FORCE)
set(WITH_OPENMP ON CACHE BOOL "OpenMP will be used" FORCE)

include(FetchContent)
FetchContent_Declare(mahi-util GIT_REPOSITORY https://github.com/mahilab/mahi-util.git)
FetchContent_MakeAvailable(mahi-util)
FetchContent_Declare(casadi GIT_REPOSITORY https://github.com/casadi/casadi.git GIT_TAG fadc86444f3c7ab824dc3f2d91d4c0cfe7f9dad5)
FetchContent_MakeAvailable(casadi)

if (MSVC)
        add_compile_options(/openmp)
endif()

set(MAHI_MPC_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/ModelGenerator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/TrajectoryGenerator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/WaypointGenerator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/ModelParameters.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/ModelControl.cpp)

file(GLOB_RECURSE MAHI_MPC_HEADERS "include/*.hpp")

add_library(mpc ${CMAKE_CURRENT_SOURCE_DIR}/src/Mahi/Mpc/ModelGenerator.cpp)
add_library(mahi::mpc ALIAS mpc)
set_target_properties(mpc PROPERTIES DEBUG_POSTFIX -d)
target_compile_features(mpc PUBLIC cxx_std_11)
set_target_properties(mpc PROPERTIES OUTPUT_NAME mpc)

target_compile_definitions(mpc PUBLIC MAHI_MPC) # for compatibility checks


install(TARGETS mpc EXPORT mpc-targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

target_sources(mpc PRIVATE ${MAHI_MPC_SRC} ${MAHI_MPC_HEADERS}) # for intellisense
target_include_directories(mpc
        PUBLIC
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(mpc casadi mahi::util)

# If this is just being included to another project, this copies casadi.dll and 
# casadi_nlpsol_ipopt.dll to the main projects build/config directory. Ideally,
# this would be abstracted for compiler (I think this will only work for MSVC), 
# but this works for now.
IF (NOT ${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_command(
        TARGET mpc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/build/_deps/casadi-build/$<CONFIG>/casadi.dll
                ${CMAKE_SOURCE_DIR}/build/$<CONFIG>/casadi.dll)       
add_custom_command(
        TARGET mpc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/build/_deps/casadi-build/$<CONFIG>/casadi_nlpsol_ipopt.dll
                ${CMAKE_SOURCE_DIR}/build/$<CONFIG>/casadi_nlpsol_ipopt.dll)       
ENDIF()

# this adds a function that can be used by executables if they are not output in the main build directory
function(link_mpc mpc_target_name)
add_custom_command(
    TARGET ${mpc_target_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/build/_deps/casadi-build/$<CONFIG>/casadi.dll
            $<TARGET_FILE_DIR:${mpc_target_name}>)     
add_custom_command(
    TARGET ${mpc_target_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/build/_deps/casadi-build/$<CONFIG>/casadi_nlpsol_ipopt.dll
            $<TARGET_FILE_DIR:${mpc_target_name}>)     
endfunction()

if(MAHI_MPC_EXAMPLES)
add_subdirectory(examples)
endif()
