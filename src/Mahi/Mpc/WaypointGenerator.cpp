//#include <MOE/Moe.hpp>
#include <Mahi/Mpc/WaypointGenerator.hpp>
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
    
    int np = m_model_parameters.waypoint_list.size();
    std::cout<<np;
    casadi::SX curr_pos_fun = SX::sym("curr_pos_fun", m_model_parameters.num_x_t);

    casadi::SX x_next = m_x + m_x_dot*m_model_parameters.step_size_t.as_seconds();
    casadi::Function F = casadi::Function("F",{m_x,m_u},{x_next},{"x","u"},{"x_next"});



//Data saving vectors
    std::vector<std::vector<double>> data;
    std::vector<double> data_line;
    std::vector<int> dof = m_model_parameters.dof;
    int dof_val;
    if(!std::count(dof.begin(),dof.end(),2)) dof_val = 3;
    else if(!std::count(dof.begin(),dof.end(),0)) dof_val = 5;
    else dof_val = 7;

//Robot state bounds
    std::vector<double> x_min, x_max;
    for(int i = 0; i<m_model_parameters.dof.size(); i++){
    x_min.push_back(m_model_parameters.x_min[dof[i]]);
    x_max.push_back(m_model_parameters.x_max[dof[i]]);
    }
    for(int i = 0; i<m_model_parameters.dof.size(); i++){
    x_min.push_back(-10*mahi::util::DEG2RAD);
    x_max.push_back(10*mahi::util::DEG2RAD);
    }
    
//Control (Alpha) initial guess and bounds 
    std::vector<double> u_init(m_model_parameters.num_u_t,0.0); 
    std::vector<double> u_min(m_model_parameters.num_u_t,0);  
    std::vector<double> u_max(m_model_parameters.num_u_t,1);  
        // Data saving information
    
//Loop for each local minima/maxima on trajectory
   for(int i = 0; i<np; i++ ){

        double time_span = m_model_parameters.waypoint_list[i+1][0] - m_model_parameters.waypoint_list[i][0];
        double num_shooting_nodes_t = time_span/m_model_parameters.step_size_t.as_seconds();
    //Initializing initial and terminal state constraints and initial guess
        std::vector<double> x0_min = m_model_parameters.waypoint_list[i];
        std::vector<double> xf_min = m_model_parameters.waypoint_list[i+1];
        
        x0_min.erase(x0_min.begin());
        xf_min.erase(xf_min.begin());

    //These likely could/ should be pointers for computation efficiency
        std::vector<double> x_init = x0_min;
        std::vector<double> x0_max = x0_min;
        std::vector<double> xf_max = xf_min;
        //std::cout<<x0_min<<std::endl;
        //std::cout<<xf_min<<std::endl;


    //Total number of variables (number of states * (time steps+1)  +  number of control inputs * time steps 
        int NV = m_model_parameters.num_x_t*(num_shooting_nodes_t+1) + m_model_parameters.num_u_t*num_shooting_nodes_t;      
    
    //Declare a variable vector to use in NLP
        casadi::MX V = casadi::MX::sym("V",NV);

    //NLP variable bounds and initial guesses
        std::vector<double> v_min,v_max,v_init;

    //Offset in V --- this is like a counter variable
        int offset=0;

    //Declare vectors for the state and control at each node, insert state and control bounds
        std::vector<casadi::MX> X, U;
        for(int k=0; k<num_shooting_nodes_t; ++k){
        //Local state
            X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x_t))));
            if(k==0){
                v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
                v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
            }
            else {
                v_min.insert(v_min.end(), x_min.begin(), x_min.end());
                v_max.insert(v_max.end(), x_max.begin(), x_max.end());
            }
            v_init.insert(v_init.end(), x_init.begin(), x_init.end());
            offset += m_model_parameters.num_x_t;
        //Local Control
            U.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_u_t))));
            v_min.insert(v_min.end(), u_min.begin(), u_min.end());
            v_max.insert(v_max.end(), u_max.begin(), u_max.end());
            v_init.insert(v_init.end(), u_init.begin(), u_init.end());
            offset += m_model_parameters.num_u_t;
        }

    //State at end
        X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x_t))));
        v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
        v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += m_model_parameters.num_x_t;

    //Objective function
        casadi::MX J = 0;
        double current_t = 0.0;

    //Constratin function and bounds
        std::vector<casadi::MX> g;
    //Loop over shooting nodes for calculating dynamics and cost function
        for(int k=0; k<num_shooting_nodes_t; ++k){
        //Temporary variable for defining the mean of activition squared
            MX temp_J = 0;
            MX averaged_J_node = 0;
        //Defining current state from xi+1 = f(xi,ui)     
            casadi::MX current_state = casadi::MX::sym("current_state", m_model_parameters.num_x_t);

            auto funcOut =\
            F({{"x",X[k]},{"u",U[k]}});
            current_state = funcOut["x_next"];
        // Save dynamic constraints
            g.push_back(current_state-X[k+1]);
        
        //Temp MX to add alpha squared for each time step
            for(int j = 0; j<m_model_parameters.num_u_t; j++){
                temp_J += (mtimes(U[k](j),U[k](j)));
            }
        //Cost function contribution of mean alpha^2 
            J +=temp_J/m_model_parameters.num_u_t;
        }

    //Vertically concatenating dynamic constraints
        casadi::MX g_vec = casadi::MX::vertcat(g);
    // NLP formulation
        casadi::MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec}};
    //Default opts for IPOPT
        casadi::Dict opts;
    //Create an NLP solver and buffers
        m_solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
   
    //Adding state and control bounds and initial guesses to problem
        std::map<std::string, DM> arg, res;
        arg["lbx"] = v_min;
        arg["ubx"] = v_max;
        arg["lbg"] = 0;
        arg["ubg"] = 0;
        arg["x0"] = v_init;
    //Solve the problem
        res = m_solver(arg);
   
    // Optimal solution of the NLP
        std::vector<double> V_opt(res.at("x"));
   
    // Data pushback for shooting nodes
        for(int k = 0; k<num_shooting_nodes_t;k++){
            data_line.clear();
            //states at node
            data_line.push_back(m_model_parameters.step_size_t.as_seconds()*k+ m_model_parameters.waypoint_list[i][0]);
            for (int i = 0; i < m_model_parameters.num_x_t; i++) data_line.push_back(V_opt[k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i]);//[) ;
            //control at node
            for (int i = 0; i < m_model_parameters.num_u_t; i++) data_line.push_back(V_opt[k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i+m_model_parameters.num_x_t]);
            data.push_back(data_line);       
        }
    std::vector<std::string> header;
    if(dof_val == 3){
    header = {"Time (s)",
                                    "EFE ref (rad)",  // "FPS ref (rad)",   "WFE ref (rad)",   "WRU ref (rad)",
                                    "EFE ref (rad/s)", //"FPS ref (rad/s)", "WFE ref (rad/s)", "WRU ref (rad/s)"
                                    "Alpha 0",
                                    "Alpha 1"};
    }
    if(dof_val == 5){
    header = {"Time (s)",
                                    "WFE ref (rad)",  // "FPS ref (rad)",   "WFE ref (rad)",   "WRU ref (rad)",
                                    "WFE ref (rad/s)", //"FPS ref (rad/s)", "WFE ref (rad/s)", "WRU ref (rad/s)"
                                    "Alpha 4",
                                    "Alpha 5"};
    }
    if(dof_val == 7){
    header = {"Time (s)",
                                    "EFE ref (rad)",  // "FPS ref (rad)",   "WFE ref (rad)",   "WRU ref (rad)",
                                    "WFE ref (rad)",
                                    "EFE ref (rad/s)",
                                    "WFE ref (rad/s)", //"FPS ref (rad/s)", "WFE ref (rad/s)", "WRU ref (rad/s)"
                                    "Alpha 0",
                                    "Alpha 1",
                                    "Alpha 4",
                                    "Alpha 5"};
    }
    mahi::util::Timestamp ts;
    std::string save_filepath = "C:/Git/fes-exo-traj-opt/deidentified_data/" + m_model_parameters.name_t + "/Opt_Traj/" + std::to_string(dof_val) + "_optimized_trajectory";
    std::cout << save_filepath << std::endl;
    mahi::util::csv_write_row(save_filepath + ".csv",header); 
    mahi::util::csv_append_rows(save_filepath + ".csv",data);   
    } // Waypoint list size - Number of IPOPT Iterations
} // Create trajectory
} // Namespace mpc
} // Namespace mahi