cmake_minimum_required(VERSION 3.13)
project(casadi_tests)

set(CMAKE_CXX_FLAGS "-openmp")

set(IPOPT_BASE_FILEPATH "C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md")
set(IPOPT_LIBRARIES "${IPOPT_BASE_FILEPATH}/lib/ipopt.dll.lib")
set(IPOPT_INCLUDE_DIRS "${IPOPT_BASE_FILEPATH}/include/coin-or")
set(Ipopt_DIR "${IPOPT_BASE_FILEPATH}/bin")
set(WITH_IPOPT ON CACHE BOOL "IPOPT will be used" FORCE)
# set(WITH_CLANG OFF CACHE BOOL "CLANG will be used" FORCE)
set(WITH_OPENMP ON CACHE BOOL "OpenMP will be used" FORCE)
#option(WITH_OPENMP "Compile with parallelization support" OFF)


include(FetchContent)
FetchContent_Declare(mahi-util GIT_REPOSITORY https://github.com/mahilab/mahi-util.git)
FetchContent_MakeAvailable(mahi-util)
FetchContent_Declare(casadi GIT_REPOSITORY https://github.com/casadi/casadi.git)
FetchContent_MakeAvailable(casadi)


if (MSVC)
        add_compile_options(/openmp)
endif()


include_directories(include)

add_executable(mainTest src/mainTest.cpp)
target_link_libraries(mainTest mahi::util casadi)

add_executable(codegen_usage src/codegen_usage.cpp)
target_link_libraries(codegen_usage mahi::util casadi)

# add_executable(nlp_codegen src/nlp_codegen.cpp)
# target_link_libraries(nlp_codegen mahi::util casadi)

add_executable(codegen_playground src/codegen_playground.cpp)
target_link_libraries(codegen_playground mahi::util casadi)

add_executable(nlp_multishoot src/nlp_multishoot.cpp)
target_link_libraries(nlp_multishoot mahi::util casadi)

add_executable(test_codegen src/test_codegen.cpp)
target_link_libraries(test_codegen mahi::util casadi)

# file(GLOB_RECURSE THEDLLS "C:/Git/casadi - clang/build/*.dll")
# foreach(currentDLL ${THEDLLS})
#     add_custom_command(TARGET nlp_multishoot POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                     ${currentDLL}
#                     "${PROJECT_SOURCE_DIR}/build/Release")
# endforeach(currentDLL)
# file(GLOB_RECURSE THEDLLS "${PROJECT_SOURCE_DIR}/build/deps/casadi-build/Release/*.dll")
# foreach(currentDLL ${THEDLLS})
#     add_custom_command(TARGET nlp_multishoot POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                     ${currentDLL}
#                     "${PROJECT_SOURCE_DIR}/build/Release")
# endforeach(currentDLL)

set(MAHI_CASADI_SRC
        src/Mahi/Casadi/ModelGenerator.cpp
        include/Mahi/Casadi/ModelGenerator.hpp)

# add_library(mahi-casadi ${MAHI_CASADI_SRC})
# target_include_directories(mahi-casadi PUBLIC include)

# file(GLOB_RECURSE MAHI_CASADI_HEADERS "include/*.hpp" "include/*.h" "include/*.inl")
# target_sources(mahi-casadi PRIVATE ${MAHI_CASADI_HEADERS}) # for intellisense
# target_include_directories(mahi-casadi
#         PUBLIC
#                 $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
#                 $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#         PRIVATE
#                 ${CMAKE_CURRENT_SOURCE_DIR}/src
# )

# target_link_libraries(mahi-casadi casadi mahi-util)

add_executable(model_generator_test src/model_generate_example.cpp ${MAHI_CASADI_SRC})
target_link_libraries(model_generator_test casadi mahi::util)

add_executable(moe_model_generator_test src/moe_model_generate_example.cpp ${MAHI_CASADI_SRC})
target_link_libraries(moe_model_generator_test casadi mahi::util)

# add_executable(model_control_test src/model_control_example.cpp ${MAHI_CASADI_SRC})
# target_link_libraries(model_control_test casadi mahi::util)