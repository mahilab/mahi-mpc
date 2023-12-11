//#include <MOE/Moe.hpp>
#include <Mahi/Mpc/TrajectoryGenerator.hpp>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;
//using namespace moe;

namespace mahi {
namespace mpc {

TrajectoryGenerator::TrajectoryGenerator(TrajectoryParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX u): 
    m_model_parameters(model_parameters),
    m_x(x),
    m_x_dot(x_dot),
    m_u(u)
{
}

TrajectoryGenerator::~TrajectoryGenerator(){

}

void TrajectoryGenerator::create_trajectory(){

    //mahi::util::print("generating a {}linear model with {} shooting nodes over {} seconds with {} states, and {} control variables", (m_model_parameters.is_linear ? "" : "non"), m_model_parameters.num_shooting_nodes, m_model_parameters.timespan.as_seconds(), m_model_parameters.num_x, m_model_parameters.num_u);
    int np = 1;
    casadi::SX curr_pos_fun = SX::sym("curr_pos_fun", m_model_parameters.num_x_t);

    // nonlinear version
    casadi::SX x_next = m_x + m_x_dot*m_model_parameters.step_size_t.as_seconds();
    casadi::Function F = casadi::Function("F",{m_x,m_u},{x_next},{"x","u"},{"x_next"});

    // total number of variables (number of states * (time steps+1)  +  number of control inputs * time steps )
    int NV = m_model_parameters.num_x_t*(m_model_parameters.num_shooting_nodes_t+1) + m_model_parameters.num_u_t*m_model_parameters.num_shooting_nodes_t;


    std::vector<double> zero_points = {0,0,0,0};

//data saving vectors
    std::vector<std::vector<double>> data;
    std::vector<double> data_line;

     //   std::array<double, 4> position_max = params.pos_limits_max_;

    std::vector<double> x_min, x_max;
    for(int i = 0; i<=m_model_parameters.dof; i++){
    x_min.push_back(m_model_parameters.x_min[i]);
    x_max.push_back(m_model_parameters.x_max[i]);
    }

    for(int i = 0; i<=m_model_parameters.dof; i++){
    x_min.push_back(-10*mahi::util::RAD2DEG);
    x_max.push_back(10*mahi::util::RAD2DEG);
    }

    std::vector<double> u_min(m_model_parameters.num_u_t,0);  
    std::vector<double> u_max(m_model_parameters.num_u_t,1);  

    for(int i = 0; i<=np; i++ ){

    // initial guess for u and x. We are guessing that all states as well as control inputs start at 0
    std::vector<double> u_init(m_model_parameters.num_u_t,0.0);
    std::vector<double> x_init(m_model_parameters.num_x_t,0.0);

    std::vector<double> x0_min(m_model_parameters.num_x_t,0.0);
    std::vector<double> x0_max(m_model_parameters.num_x_t,0.0);

    // we don't set a found on final condition
    std::vector<double> xf_min(m_model_parameters.num_x_t,-casadi::inf);
    std::vector<double> xf_max(m_model_parameters.num_x_t,casadi::inf);

    // inital condition for state. This is essentially saying that we start at state of 0 for both position and velocity.
    //Loop here - need to pull in the csv that has randomized setpoints.

    //if(i == 0){
    xf_min = {-70.0 * util::DEG2RAD,0};// -70* util::DEG2RAD, 0, 0};
    xf_max = {-70.0 *util::DEG2RAD, 0};// 70* util::DEG2RAD, 0, 0}; // desired_waypoints[i+1];
  //  }
    // else{
    // x0_min = desired_waypoints[i];
    // x0_max = desired_waypoints[i];

    // xf_min = desired_waypoints[i+1];
    // xf_max = desired_waypoints[i+1];
    // }

    // declare a variable vector to use in NLP
    casadi::MX V = casadi::MX::sym("V",NV);

    // NLP variable bounds and initial guesses
    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    // declare vectors for the state and control at each node
    std::vector<casadi::MX> X, U;
    for(int k=0; k<m_model_parameters.num_shooting_nodes_t; ++k){
        // Local state
        X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x_t))));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += m_model_parameters.num_x_t;

        // Local Control
        U.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_u_t))));
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
        offset += m_model_parameters.num_u_t;
    }

    // State at end
    X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x_t))));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += m_model_parameters.num_x_t;

    // Objective function
    casadi::MX J = 0;
    double current_t = 0.0;
    // Constratin function and bounds
    std::vector<casadi::MX> g;

    int traj_size = m_model_parameters.num_shooting_nodes_t*m_model_parameters.num_x_t;
//WHAT IS THIS FOR?
//traj_size += m_model_parameters.num_u;            // u_init

    // create a vector of symbolic variables that we will import. This is of size (num_time_steps * num_states)
    casadi::MX traj = casadi::MX::sym("traj", traj_size); // traj_size
    //DO I STILL NEED THIS? What does this part do? 

    //     int start_u_init     = end_Rm;
    //     int end_u_init       = start_u_init + m_model_parameters.num_u;
    //     u_init_in  = reshape(traj(casadi::Slice(start_u_init,end_u_init)),m_model_parameters.num_u,1);
    // }


    // Loop over shooting nodes
    for(int k=0; k<m_model_parameters.num_shooting_nodes_t; ++k){
        // iterate through dynamics using either linearized or nonlinear dynamics
        casadi::MX current_state = casadi::MX::sym("current_state", m_model_parameters.num_x_t); // traj_size

        auto funcOut =\
        F({{"x",X[k]},{"u",U[k]}});
        current_state = funcOut["x_next"];
        
        // Save continuity constraints
        g.push_back(current_state-X[k+1]);

        // Add objective function contribution on change of input
        //J += mtimes(U[k],U[k]);
        MX temp_J, averaged_J_node;
        for(int j = 0; j<=m_model_parameters.num_u_t; j++){
             temp_J += mtimes(U[k](j),U[k](j));
        }
        averaged_J_node = temp_J/m_model_parameters.num_u_t;
        J += averaged_J_node;
    }

    // NLP
    casadi::MX g_vec = casadi::MX::vertcat(g);
    casadi::MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec}, {"p", traj}};

    // leave default but can add later at runtime
    casadi::Dict opts;

    // Create an NLP solver and buffers
    m_solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
    std::map<std::string, DM> arg, res;
    // // Bounds and initial guess
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;
    //Solve the problem
    res = m_solver(arg);
    // Optimal solution of the NLP
    std::vector<double> V_opt(res.at("x"));
    std::cout<<V_opt<<std::endl;
//     for(int k = 0; k<=m_model_parameters.np_; k++){
//             //states at node
//             for (int i = 0; i < m_model_parameters.num_x_t; i++) data_line.push_back(V_opt[k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i]);//[) ;
//             //control at node
//             for (int i = 0; i < m_model_parameters.num_x_t; i++) data_line.push_back(V_opt[k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i+m_model_parameters.num_x_t]);
//             data.push_back(data_line);
//     }
//    //
//       // data saving information
//     std::vector<std::string> header = {"Time (s)", 
//                                     "EFE ref (rad)",   "FPS ref (rad)",   "WFE ref (rad)",   "WRU ref (rad)",
//                                     "EFE ref (rad/s)", "FPS ref (rad/s)", "WFE ref (rad/s)", "WRU ref (rad/s)"};
//     mahi::util::Timestamp ts;
//     std::string save_filepath = "C:/Git/fes-exo-traj-opt/deidentified_data/S" + m_model_parameters.name_t + "/optimized_trajectory" + ts.yyyy_mm_dd_hh_mm_ss();
//     std::cout << save_filepath << std::endl;

//     mahi::util::csv_write_row(save_filepath + ".csv",header);
//     mahi::util::csv_append_rows(save_filepath + ".csv",data);

}

// std::vector<std::vector<double>> get_waypoint_list(std::string traj_filepath, int n_rows, std::vector<int> dof, std::vector<double> zero_positions){
//     std::vector<std::vector<double>> file_data(n_rows, std::vector<double>(dof.size()));
//     mahi::util::csv_read_rows(traj_filepath, file_data, 1, 0); 

//     //std::vector<mahi::robo::WayPoint> traj_waypoints;
//     for (const auto &waypoint : file_data){
//         std::vector<std::vector<double>> waypoint_list;
//         std::vector<double> positions(waypoint.begin()+1,waypoint.end());
//         for (size_t i = 0; i < positions.size(); i++){
//             positions[i]*= mahi::util::DEG2RAD;    
//             if (!std::count(dof.begin(),dof.end(),i)){
//                 positions[i] = zero_positions[i];
//             }
//         waypoint_list.push_back(positions);
//         } 
//     return waypoint_list;
//     }
    
// }

// void ModelGenerator::generate_c_code(){
//     m_c_file_filepath = m_model_parameters.name + ".c";
//     std::cout << "generating c code to filepath " << m_c_file_filepath << std::endl;
//     m_solver.generate_dependencies(m_c_file_filepath);
// }


// void ModelGenerator::compile_model(){
//     m_model_parameters.dll_filepath = m_model_parameters.name + ".so";
//     std::cout << "generating dll at filepath " << m_model_parameters.dll_filepath << std::endl;
//     int flag = system(("gcc -fPIC -shared -O1 " + m_c_file_filepath + " -o " + m_model_parameters.dll_filepath).c_str());
//     casadi_assert(flag==0, "Compilation failed");
//     save_param_file();
// }

// void ModelGenerator::save_param_file(){
//     nlohmann::json j;

//     j["model"] = m_model_parameters;

//     std::ofstream file1(m_model_parameters.name + ".json");
//     if (file1.is_open())
//         file1 << j;
//     file1.close();
}
} // namespace mpc
} // namespace mahi