#pragma once

#include <string>
#include <vector>
#include <Mahi/Casadi/ModelParameters.hpp>
#include <Mahi/Util/Timing/Time.hpp>
#include <casadi/casadi.hpp>
#include <mutex>

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

    ModelControl(std::string model_name, std::vector<double> Q= {}, std::vector<double> R = {}, std::vector<double> Rm = {}, casadi::Dict solver_opts = casadi::Dict());
    ~ModelControl();

    ModelParameters model_parameters;
    std::vector<ControlResult> control_results;

    void calc_u(mahi::util::Time time,const std::vector<double>& state, const std::vector<double>& control,std::vector<double> traj);
    void load_model(const std::string& model_name);

    ControlResult control_at_time(mahi::util::Time time);

    void start_calc();
    void stop_calc();

    void set_state(mahi::util::Time time,const std::vector<double>& state, const std::vector<double>& control, std::vector<double> traj);

    void update_weights(std::vector<double> Q = {},std::vector<double> R = {},std::vector<double> Rm = {});

    void update_control_limits(std::vector<double> u_min, std::vector<double> u_max);
private:
    casadi::Dict m_solver_opts;
    mahi::util::Time curr_time;
    bool m_is_linear;
    casadi::Function m_solver;
    std::map<std::string, casadi::DM> m_solver_args;
    std::map<std::string, casadi::DM> m_solver_result;

    std::vector<double> m_Q;
    std::vector<double> m_R;
    std::vector<double> m_Rm;

    casadi::Function get_A;
    casadi::Function get_B;
    casadi::Function get_x_dot;

    std::vector<double> v_min;
    std::vector<double> v_max;

    std::atomic<bool> m_stop = false;

    mahi::util::Time m_time;
    std::vector<double> m_state;
    std::vector<double> m_control;
    std::vector<double> m_traj;
    
    std::mutex m_state_mutex;
    std::mutex m_output_mutex;
    std::mutex m_control_limits_mutex;

    std::atomic<bool> m_done_calcing = true;

    void format_outputs(std::vector<double> opt_output);
    
    int m_num_control_inputs_saved = 0;
};