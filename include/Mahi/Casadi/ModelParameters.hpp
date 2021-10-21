#pragma once

#include <Mahi/Util/Timing/Time.hpp>
#include <Mahi/Util.hpp>
#include <vector>
#include <string>

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

void to_json(mahi::util::json& j, const ModelParameters& p);

void from_json(const mahi::util::json& j, ModelParameters& p);