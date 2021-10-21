#pragma once

#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Casadi/ModelParameters.hpp>

class ModelGenerator
{
private:
    ModelParameters m_model_parameters;
    casadi::SX m_x; // number of control inputs
    casadi::SX m_x_dot; // number of control inputs
    casadi::SX m_u; // number of control inputs

    std::string m_c_file_filepath;
    casadi::Function m_solver;
    
public:
    
    ModelGenerator(ModelParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX u);
    ~ModelGenerator();
    void create_model();
    void generate_c_code();
    void generate_linear_functions(casadi::Function get_A, casadi::Function get_B, casadi::Function get_x_dot_init);
    void compile_model();
    void save_param_file();

};