
# add_executable(model_generator_test model_generate_example.cpp)
# target_link_libraries(model_generator_test mahi::mpc)

add_executable(model_generate ex_model_generate.cpp)
target_link_libraries(model_generate mahi::mpc)
link_mpc(model_generate)

# add_executable(model_control ex_model_control.cpp)
# target_link_libraries(model_control mahi::mpc)

add_executable(model_control model_control_example.cpp)
target_link_libraries(model_control mahi::mpc)
link_mpc(model_generate)

add_executable(thread_model_control thread_model_control_example.cpp)
target_link_libraries(thread_model_control mahi::mpc)
link_mpc(model_generate)

# add_executable(moe_model_generator_test moe_model_generate_example.cpp ${MAHI_CASADI_SRC})
# target_link_libraries(moe_model_generator_test casadi mahi::util)