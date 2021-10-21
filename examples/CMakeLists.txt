
add_executable(mainTest mainTest.cpp)
target_link_libraries(mainTest mahi::util casadi)

add_executable(model_generator_test model_generate_example.cpp)
target_link_libraries(model_generator_test mahi::mpc)

add_executable(model_generate ex_model_generate.cpp)
target_link_libraries(model_generate mahi::mpc)

add_executable(model_control model_control_example.cpp)
target_link_libraries(model_control mahi::mpc)


# add_executable(moe_model_generator_test moe_model_generate_example.cpp ${MAHI_CASADI_SRC})
# target_link_libraries(moe_model_generator_test casadi mahi::util)