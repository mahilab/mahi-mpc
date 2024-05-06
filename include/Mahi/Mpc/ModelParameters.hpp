#pragma once

#include <Mahi/Util/Timing/Time.hpp>
#include <Mahi/Util.hpp>
#include <vector>
#include <string>

namespace mahi {
namespace mpc {

struct ModelParameters{
    ModelParameters(std::string name_, int num_x_, int num_u_, mahi::util::Time step_size_, size_t num_shooting_nodes_, bool is_linear_, std::vector<double> u_min_ = {}, std::vector<double> u_max_ = {}, std::vector<double> x_min_ = {},  std::vector<double> x_max_ = {});

    ModelParameters() {} ;

    std::string name;           // name of the model for useful outputs
    mahi::util::Time timespan;  // duration of the MPC time period.
    mahi::util::Time step_size; // step size for MPC calculations 
    int num_x;               // number of states
    int num_u;               // number of control inputs
    int num_shooting_nodes;  // number of shooting_ndes
    std::vector<double> x_min;  // minimum of states
    std::vector<double> u_min;  // minimum of control inputs
    std::vector<double> x_max;  // maximum of states
    std::vector<double> u_max;  // maximum of control inputs
    std::string dll_filepath;   // filepath to dll
    bool is_linear;             // whether or not this is a linearized model
};

struct TrajectoryParameters{
    TrajectoryParameters(std::string name_t_, int num_x_t_, int num_u_t_, mahi::util::Time step_size_t_, size_t num_shooting_nodes_, std::vector<int> dof, int np, std::array<double, 4Ui64> x_min_, std::array<double, 4Ui64> x_max_, std::vector<std::vector<double>> waypoint_list_, std::vector<bool> muscles_enabled_);
    TrajectoryParameters() {} ;

    std::string name_t;           // name of the model for useful outputs
    mahi::util::Time timespan_t;  // duration of the MPC time period.
    mahi::util::Time step_size_t; // step size for MPC calculations 
    int num_x_t;               // number of states
    int num_u_t;               // number of control inputs
    int num_shooting_nodes_t;  // number of shooting_ndes
    std::vector<int> dof;
    int np;
    std::array<double, 4Ui64> x_min;
    std::array<double, 4Ui64> x_max;
    std::vector<std::vector<double>> waypoint_list;
    std::vector<bool> muscles_enabled;
};

void to_json(mahi::util::json& j, const ModelParameters& p);

void from_json(const mahi::util::json& j, ModelParameters& p);

} // namespace mpc
} // namespace mahi