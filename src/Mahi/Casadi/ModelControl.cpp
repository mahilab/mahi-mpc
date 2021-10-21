#include <Mahi/Util.hpp>
#include <Mahi/Casadi/ModelControl.hpp>

ModelControl::ModelControl(std::string model_name, casadi::Dict solver_opts):
m_solver_opts(solver_opts)
{
    load_model(model_name);
}

void ModelControl::load_model(const std::string& model_name){

    mahi::util::json j;
    std::ifstream file(model_name + ".json");
    file >> j;
    model_parameters = j["model"].get<ModelParameters>();

    // std::cout << "loaded model\n";

    // default to zero, but will change this after one iteration
    std::vector<double> x_init(model_parameters.num_x,0);
    std::vector<double> u_init(model_parameters.num_u,0);

    // std::vector<double> v_min;
    // std::vector<double> v_max;
    std::vector<double> v_init;

    //declare vectors for the state and control at each node
    for(int k=0; k<model_parameters.num_shooting_nodes; ++k){
        v_min.insert(v_min.end(), model_parameters.x_min.begin(), model_parameters.x_min.end());
        v_max.insert(v_max.end(), model_parameters.x_max.begin(), model_parameters.x_max.end());
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());

        // Local Control
        v_min.insert(v_min.end(), model_parameters.u_min.begin(), model_parameters.u_min.end());
        v_max.insert(v_max.end(), model_parameters.u_max.begin(), model_parameters.u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    }

    v_min.insert(v_min.end(), model_parameters.x_min.begin(), model_parameters.x_min.end());
    v_max.insert(v_max.end(), model_parameters.x_max.begin(), model_parameters.x_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    // Bounds and initial guess
    casadi::Dict solver_opts;

    solver_opts["ipopt.tol"] = 1e-5;
    solver_opts["ipopt.max_iter"] = 200;
    solver_opts["ipopt.linear_solver"] = "mumps"; // ma27
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["print_time"] = 0;
    solver_opts["ipopt.sb"] = "yes";

    std::string dll = model_parameters.dll_filepath;
    m_solver = nlpsol("nlpsol", "ipopt", dll, solver_opts);

    m_solver_args["lbx"] = v_min;
    m_solver_args["ubx"] = v_max;
    m_solver_args["lbg"] = 0;
    m_solver_args["ubg"] = 0;
    m_solver_args["x0"] = v_init;

    get_A = casadi::external(model_parameters.name + "_get_A",model_parameters.name + "_linear_functions.so");
    get_B = casadi::external(model_parameters.name + "_get_B",model_parameters.name + "_linear_functions.so");
    get_x_dot = casadi::external(model_parameters.name + "_get_x_dot_init",model_parameters.name + "_linear_functions.so");
}

std::vector<ModelControl::ControlResult> ModelControl::calc_u(mahi::util::Time control_time,const std::vector<double>& state, const std::vector<double>& control, std::vector<double> traj){
    curr_time = control_time;

    
    // set A and B if necessary
    if (model_parameters.is_linear){
        std::vector<casadi::DM> linear_args = {state, control};
        std::vector<double> A_res(get_A(linear_args)[0]);
        std::vector<double> B_res(get_B(linear_args)[0]);
        std::vector<double> x_dot_init_res(get_x_dot(linear_args)[0]);

        for (const auto &a_ : A_res) traj.push_back(a_);
        for (const auto &b_ : B_res) traj.push_back(b_);
        for (const auto &x_dot_init_ : x_dot_init_res) traj.push_back(x_dot_init_);
        for (const auto &x_next_ : state) traj.push_back(x_next_);
        for (const auto &u_ : control) traj.push_back(u_);
    }

    // set input args
    m_solver_args["p"]  = traj;

    std::copy(state.begin(),state.end(),v_min.begin());
    std::copy(state.begin(),state.end(),v_max.begin());

    m_solver_args["lbx"] = v_min;
    m_solver_args["ubx"] = v_max;

    // solve
    m_solver_result = m_solver(m_solver_args);
    std::vector<double> result_vec(m_solver_result.at("x"));

    m_solver_args["x0"] = result_vec;

    format_outputs(result_vec);
    // for (const auto &i : control_results){
    //     std::cout << i.time << ": " << i.x_est << ", " << i.u << std::endl;
    // }
    
    return control_results;
}

std::vector<ModelControl::ControlResult> ModelControl::format_outputs(std::vector<double> opt_output){
    control_results.clear();
    auto nx = model_parameters.num_x;
    auto nu = model_parameters.num_u;
    for (size_t i = 0; i < model_parameters.num_shooting_nodes; i++){
        std::vector<double> new_x_vec(opt_output.begin() + i*(nx+nu),      opt_output.begin() + i*(nx+nu) + nx);
        std::vector<double> new_u_vec(opt_output.begin() + i*(nx+nu) + nx, opt_output.begin() + i*(nx+nu) + nx + nu);
        control_results.emplace_back(mahi::util::seconds(curr_time.as_seconds() + model_parameters.step_size.as_seconds()*i),
                                     new_x_vec,
                                     new_u_vec);
    }
    return control_results;
}

ModelControl::ControlResult ModelControl::control_at_time(mahi::util::Time time){
    int i = 0;
    while (control_results[i].time < time) i++;
    return control_results[i];
}