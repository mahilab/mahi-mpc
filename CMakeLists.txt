cmake_minimum_required(VERSION 3.13)
project(casadi_tests)

include(FetchContent)
FetchContent_Declare(mahi-util GIT_REPOSITORY https://github.com/mahilab/mahi-util.git)
FetchContent_MakeAvailable(mahi-util)
FetchContent_Declare(casadi GIT_REPOSITORY https://github.com/casadi/casadi.git)
FetchContent_MakeAvailable(casadi)
# set(IPOPT_INCLUDEDIR "C:/Users/jeff/AppData/Local/Temp/Ipopt-3.14.2-win64-msvs2019-md.zip/Ipopt-3.14.2-win64-msvs2019-md/include/coin-or")
# set(IPOPT_LIBRARIES "C:/Users/jeff/AppData/Local/Temp/Ipopt-3.14.2-win64-msvs2019-md.zip/Ipopt-3.14.2-win64-msvs2019-md/lib" )
# set(IPOPT_LIBRARY_DIRS "C:/Users/jeff/AppData/Local/Temp/Ipopt-3.14.2-win64-msvs2019-md.zip/Ipopt-3.14.2-win64-msvs2019-md")

# set(WITH_IPOPT OFF CACHE INTERNAL "")
# set(WITH_LAPACK OFF CACHE INTERNAL "")
# set(WITH_QPOASES OFF CACHE INTERNAL "")
# set(WITH_COMMON OFF CACHE INTERNAL "")
# FetchContent_Declare(
#         casadi 
#         GIT_REPOSITORY https://github.com/casadi/casadi.git)
#         #CMAKE_PREFIX_PATH="C:/"
#         #CMAKE_ARGS -DWITH_IPOPT=OFF)
# FetchContent_MakeAvailable(casadi)




include_directories(include)

#include_directories(C:/Git/casadi)
#include_directories(C:/MatlabHelpers/Casadi/casadi-windows-matlabR2016a-v3.5.5/include/casadi)
#include_directories(C:/MatlabHelpers/Casadi/casadi-windows-matlabR2016a-v3.5.5/include/casadi/core)
#include_directories(C:/MatlabHelpers/Casadi/casadi-windows-matlabR2016a-v3.5.5/include/casadi/core/runtime)

add_executable(mainTest src/mainTest.cpp)
target_link_libraries(mainTest mahi::util casadi)

add_executable(codegen_usage src/codegen_usage.cpp)
target_link_libraries(codegen_usage mahi::util casadi)

add_executable(nlp_codegen src/nlp_codegen.cpp)
target_link_libraries(nlp_codegen mahi::util casadi)

add_executable(codegen_playground src/codegen_playground.cpp)
target_link_libraries(codegen_playground mahi::util casadi)

add_executable(nlp_multishoot src/nlp_multishoot.cpp)
target_link_libraries(nlp_multishoot mahi::util casadi)

add_executable(test_codegen src/test_codegen.cpp)
target_link_libraries(test_codegen mahi::util casadi)

# file(GLOB_RECURSE THEDLLS "C:/Git/casadi/build/Release/*.dll")
# foreach(currentDLL ${THEDLLS})
#     add_custom_command(TARGET mainTest POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                     ${currentDLL}
#                     "${PROJECT_SOURCE_DIR}/build/Release")
# endforeach(currentDLL)

file(GLOB_RECURSE THEDLLS "C:/Git/casadi/build/Release/*.dll")
foreach(currentDLL ${THEDLLS})
    add_custom_command(TARGET nlp_multishoot POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${currentDLL}
                    "${PROJECT_SOURCE_DIR}/build/Release")
endforeach(currentDLL)

# add_library(multishoot_library SHARED "build/Release/multishoot_c.c")


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

# add_executable(model_control_test src/model_control_example.cpp ${MAHI_CASADI_SRC})
# target_link_libraries(model_control_test casadi mahi::util)