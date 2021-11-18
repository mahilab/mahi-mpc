#include <Mahi/Util.hpp>
#include <Mahi/Casadi/ModelControl.hpp>


ModelControl::ModelControl(std::string model_name, std::vector<double> Q, std::vector<double> R, casadi::Dict solver_opts):
m_Q(Q),
m_R(R),
m_solver_opts(solver_opts)
{
    std::cout << m_R << std::endl;
    load_model(model_name);
}

ModelControl::~ModelControl(){
    m_stop = true;
    while(!m_done_calcing){}
}

void ModelControl::load_model(const std::string& model_name){

    mahi::util::json j;
    std::ifstream file(model_name + ".json");
    file >> j;
    model_parameters = j["model"].get<ModelParameters>();

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

void ModelControl::set_state(mahi::util::Time time,const std::vector<double>& state, const std::vector<double>& control, std::vector<double> traj){
    std::lock_guard<std::mutex> lg(m_state_mutex);
    m_time = time;
    m_state = state;
    m_control = control;
    m_traj = traj;
}

void ModelControl::start_calc(){
    m_stop = false;
    m_done_calcing = false;
    
    std::thread new_thread([this]{
        mahi::util::Time local_time;
        std::vector<double> local_state;
        std::vector<double> local_control;
        std::vector<double> local_traj;
        
        std::vector<double> sim_times;
        std::vector<double> calc_times;
        while(!m_stop){
            mahi::util::Clock timing_clock;
            {
                std::lock_guard<std::mutex> lg(m_state_mutex);
                local_time = m_time;
                local_state = m_state;
                local_control = m_control;
                local_traj = m_traj;
            }
            calc_u(local_time, local_state, local_control, local_traj);
            sim_times.push_back(local_time.as_milliseconds());
            calc_times.push_back(timing_clock.get_elapsed_time().as_milliseconds());
        }
        for (size_t i = 0; i < sim_times.size(); i++){
            mahi::util::print("sim time: {}, calc_times: {}", sim_times[i],calc_times[i]);
        }
        m_done_calcing = true;
    });
    new_thread.detach();
}

void ModelControl::stop_calc(){m_stop = true;}

void ModelControl::calc_u(mahi::util::Time control_time,const std::vector<double>& state, const std::vector<double>& control, std::vector<double> traj){

    curr_time = control_time;

    // static std::vector<double> Q = {100, 50, 50, 50, 0.1, 0.1, 0.1, 0.1};
	// static std::vector<double> R = {0.1, 1.0, 1.0, 1.0};

    // static std::vector<double> Q = {10, 1, 1, 10, 5, 5, 5, 5};
    // static std::vector<double> R = {0,0,0,0};

    for (const auto &q : m_Q) traj.push_back(q);
    for (const auto &r : m_R) traj.push_back(r);
    
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

    auto nx = model_parameters.num_x;
    auto nu = model_parameters.num_u;

    std::copy(state.begin(),state.end(),v_min.begin());
    std::copy(state.begin(),state.end(),v_max.begin());    

    m_solver_args["lbx"] = v_min;
    m_solver_args["ubx"] = v_max;
    // solve
    m_solver_result = m_solver(m_solver_args);
    std::vector<double> result_vec(m_solver_result.at("x"));

    m_solver_args["x0"] = result_vec;

    format_outputs(result_vec);

    {
        std::lock_guard<std::mutex> lg(m_output_mutex);
        for (size_t i = 0; i < m_num_control_inputs_saved; i++){
            std::copy(control_results[i].u.begin(),control_results[i].u.end(),v_min.begin()+ i*(nx+nu) + nx);
            std::copy(control_results[i].u.begin(),control_results[i].u.end(),v_max.begin()+ i*(nx+nu) + nx);
        }
    }
}

void ModelControl::format_outputs(std::vector<double> opt_output){
    std::vector<ModelControl::ControlResult> control_results_local;
    auto nx = model_parameters.num_x;
    auto nu = model_parameters.num_u;
    for (size_t i = 0; i < model_parameters.num_shooting_nodes; i++){
        std::vector<double> new_x_vec(opt_output.begin() + i*(nx+nu),      opt_output.begin() + i*(nx+nu) + nx);
        std::vector<double> new_u_vec(opt_output.begin() + i*(nx+nu) + nx, opt_output.begin() + i*(nx+nu) + nx + nu);
        control_results_local.emplace_back(mahi::util::seconds(curr_time.as_seconds() + model_parameters.step_size.as_seconds()*i),
                                     new_x_vec,
                                     new_u_vec);
    }
    
    {
        std::lock_guard<std::mutex> lg(m_output_mutex);
        control_results = control_results_local;
    }
}

ModelControl::ControlResult ModelControl::control_at_time(mahi::util::Time time){
    std::lock_guard<std::mutex> lg(m_output_mutex);
    int i = 0;
    while (control_results[i].time < time && i < control_results.size()) i++;
    return (i == 0) ? control_results[0] : control_results[i-1];
}