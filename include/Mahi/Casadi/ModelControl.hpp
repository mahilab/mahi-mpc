#pragma once

#include <string>
#include <vector>
#include <Mahi/Casadi/ModelParameters.hpp>
#include <Mahi/Util/Timing/Time.hpp>
#include <casadi/casadi.hpp>

class ModelControl
{
public:
    struct ControlResult{
        ControlResult(mahi::util::Time time_,std::vector<double> x_est_,std::vector<double> u_):
            time(time_),
            x_est(x_est_),
            u(u_){}
        mahi::util::Time time;
        std::vector<double> x_est;
        std::vector<double> u;
    };

    ModelControl(std::string model_name, casadi::Dict solver_opts);

    ModelParameters model_parameters;
    std::vector<ControlResult> control_results;

    std::vector<ControlResult> calc_u(mahi::util::Time control_time,const std::vector<double>& state, const std::vector<double>& control,std::vector<double> traj);
    void load_model(const std::string& model_name);

    ControlResult control_at_time(mahi::util::Time time);
private:
    casadi::Dict m_solver_opts;
    mahi::util::Time curr_time;
    bool m_is_linear;
    casadi::Function m_solver;
    std::map<std::string, casadi::DM> m_solver_args;
    std::map<std::string, casadi::DM> m_solver_result;

    casadi::Function get_A;
    casadi::Function get_B;
    casadi::Function get_x_dot;

    std::vector<double> v_min;
    std::vector<double> v_max;

    std::vector<ControlResult> format_outputs(std::vector<double> opt_output);
    

};