#pragma once

#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Mpc/ModelParameters.hpp>

namespace mahi {
namespace mpc {

class TrajectoryGenerator
{
private:
    TrajectoryParameters m_model_parameters;
    casadi::SX m_x; // number of state inputs
    casadi::SX m_x_dot; // number of state inputs
    casadi::SX m_u; // number of control inputs
    //std::string m_c_file_filepath;
    casadi::Function m_solver;
    
public:
    
    TrajectoryGenerator(TrajectoryParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX u);
    ~TrajectoryGenerator();
    void create_trajectory();
    // void generate_c_code();
   
    // void compile_model();
    // void save_param_file();

};

} // namespace mpc
} // namespace mahi