#pragma once

#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Mpc/ModelParameters.hpp>

namespace mahi {
namespace mpc {

class TrajectoryGenerator
{
private:
    ModelParameters m_model_parameters;
    casadi::SX m_x; // number of state inputs
    casadi::SX m_x_dot; // number of state inputs
    casadi::SX m_u; // number of control inputs
    int 
    //std::string m_c_file_filepath;
    casadi::Function m_solver;
    
public:
    
    TrajectoryGenerator(ModelParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX , int np);
    ~TrajectoryGenerator();
    // void create_model();
    // void generate_c_code();
   
    // void compile_model();
    // void save_param_file();

};

} // namespace mpc
} // namespace mahi