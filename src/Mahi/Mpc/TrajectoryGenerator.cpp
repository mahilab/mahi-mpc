//#include <MOE/Moe.hpp>
#include <Mahi/Mpc/TrajectoryGenerator.hpp>
#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>
#include <string>

using namespace casadi;
using namespace std;
//using namespace moe;

using namespace mahi;
using namespace mpc;

TrajectoryGenerator::TrajectoryGenerator(TrajectoryParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX u, std::vector<int> cf_params, int trial_num): 
    m_model_parameters(model_parameters),
    m_x(x),
    m_x_dot(x_dot),
    m_u(u),
    m_cf_params(cf_params),
    m_trial_num(trial_num)
{
}

TrajectoryGenerator::~TrajectoryGenerator(){
}

void TrajectoryGenerator::create_trajectory(){

    int np =  m_model_parameters.waypoint_list.size() - 1;
    casadi::SX curr_pos_fun = SX::sym("curr_pos_fun", m_model_parameters.num_x_t);
    casadi::SX des_pos_fun = SX::sym("des_pos_fun", m_model_parameters.num_x_t);
    casadi::SX err_out = curr_pos_fun - des_pos_fun; 
    casadi::Function err_fun = casadi::Function("err_fun",{curr_pos_fun, des_pos_fun},{err_out},{"curr_pos_fun", "des_pos_fun"},{"err_out"});
    casadi::SX x_next = m_x + m_x_dot*m_model_parameters.step_size_t.as_seconds();
    casadi::Function F = casadi::Function("F",{m_x,m_u},{x_next},{"x","u"},{"x_next"});

//Data saving vectors
    std::vector<std::vector<double>> data;
    int trial_num = m_trial_num;
    std::vector<double> data_line;
    std::vector<int> dof = m_model_parameters.dof;
   
    std::string dof_string = "";
    std::string cf_string = "";
    for (auto d : dof) dof_string += std::to_string(d);
    for (auto d : m_cf_params) cf_string += std::to_string(d);

//Robot state bounds
    std::vector<double> x_min, x_max;
    for(int i = 0; i<m_model_parameters.dof.size(); i++){
    x_min.push_back(m_model_parameters.x_min[dof[i]]);
    x_max.push_back(m_model_parameters.x_max[dof[i]]);
    }
    for(int i = 0; i<m_model_parameters.dof.size(); i++){
    x_min.push_back(-45*mahi::util::DEG2RAD);
    x_max.push_back(45*mahi::util::DEG2RAD);
    }
    
//Control (Alpha) initial guess and bounds 
    std::vector<double> u_init(m_model_parameters.num_u_t,0.0); 
    std::vector<double> u_min(m_model_parameters.num_u_t,0.0);  
    std::vector<double> u_max(m_model_parameters.num_u_t,1.0);  

// Data saving information
    std::vector<double> Q_in, R_in, Rm_in;
    std::vector<double> Q_active, R_active, Rm_active;
    std::vector<double> Q_, R_, Rm_;

    mahi::util::json mpc_params;

    std::string mpc_params_filepath = "C:/Git/fes-exo-traj-opt/trajectories/gains.csv";


    std::vector<std::vector<double>> mpc_file_data(1, std::vector<double>(24));
    mahi::util::csv_read_rows(mpc_params_filepath, mpc_file_data, 1+ trial_num, 0); 
    for (auto &waypoint : mpc_file_data){
        for(int i =0; i<=7;i++){
            Q_in.push_back(waypoint[i]);
            R_in.push_back(waypoint[i+8]);
            Rm_in.push_back(waypoint[i+16]);
        }
    }  

    for(int i = 0; i<dof.size(); i++) Q_active.push_back(Q_in[dof[i]]);
    for(int i = 0; i<dof.size(); i++) Q_active.push_back(Q_in[4+dof[i]]);
    for(int i = 0; i <m_model_parameters.muscles_enabled.size(); i++){ 
        if(m_model_parameters.muscles_enabled[i]){
                R_active.push_back(R_in[i]);
        }
    }

    for(int i = 0; i <m_model_parameters.muscles_enabled.size(); i++){ 
        if(m_model_parameters.muscles_enabled[i]){
                Rm_active.push_back(Rm_in[i]);
        }
    }
   
    for(int i = 0; i<Q_active.size(); i++){
        for(int j = 0; j<Q_active.size(); j++){
            if(i==j){
                Q_.push_back(Q_active[i]);
            }
            else{
                Q_.push_back(0);
            }
        }
    }
    for(int i = 0; i<R_active.size(); i++){
        for(int j = 0; j<R_active.size(); j++){
            if(i==j){
                R_.push_back(R_active[i]);
                Rm_.push_back(Rm_active[i]);
            }
            else{
                R_.push_back(0);
                Rm_.push_back(0);
            }
        }
    }
 
// create a vector of symbolic variables that we will import. This includes x_init, desired cost function weights(Q,R,Rm if desired), and u_init  
    int params_size_sym = 0;
    params_size_sym += m_model_parameters.num_x_t;
    if (std::count(m_cf_params.begin(),m_cf_params.end(),0)) params_size_sym += m_model_parameters.num_x_t*m_model_parameters.num_x_t;  // Q
    if (std::count(m_cf_params.begin(),m_cf_params.end(),1)) params_size_sym += m_model_parameters.num_u_t*m_model_parameters.num_u_t;  // R
    if (std::count(m_cf_params.begin(),m_cf_params.end(),2)) params_size_sym += m_model_parameters.num_u_t*m_model_parameters.num_u_t; // Rm
    params_size_sym += m_model_parameters.num_u_t;  

//Determine where symbolic parameters lie in vector and slice 
    int start_Q =  (int)(m_model_parameters.num_x_t);  
    int start_R = start_Q; int start_Rm = start_Q; int end_Rm = start_Q;

    if (std::count(m_cf_params.begin(),m_cf_params.end(),0)){
        start_R += m_model_parameters.num_x_t*m_model_parameters.num_x_t;
        start_Rm = start_R;
        end_Rm = start_Rm;
    }
    if (std::count(m_cf_params.begin(),m_cf_params.end(),1)){
        start_Rm += m_model_parameters.num_x_t*m_model_parameters.num_x_t;
        end_Rm = start_Rm;
    }
    if (std::count(m_cf_params.begin(),m_cf_params.end(),2)){
        end_Rm   +=  m_model_parameters.num_u_t*m_model_parameters.num_u_t;
    }
    casadi::MX params_sym = casadi::MX::sym("params_sym",params_size_sym);
    casadi::MX Q, R, Rm; 
    if (std::count(m_cf_params.begin(),m_cf_params.end(),0)) Q  = reshape(params_sym(casadi::Slice(start_Q,start_R)),m_model_parameters.num_x_t,m_model_parameters.num_x_t);
    if (std::count(m_cf_params.begin(),m_cf_params.end(),1)) R  = reshape(params_sym(casadi::Slice(start_R,start_Rm)),m_model_parameters.num_u_t,m_model_parameters.num_u_t);
    if (std::count(m_cf_params.begin(),m_cf_params.end(),2)) Rm = reshape(params_sym(casadi::Slice(start_Rm,end_Rm)),m_model_parameters.num_u_t,m_model_parameters.num_u_t);

    int start_u_init     = end_Rm;
    int end_u_init       = start_u_init + m_model_parameters.num_u_t;
    casadi::MX u_init_in  = reshape(params_sym(casadi::Slice(start_u_init,end_u_init)),m_model_parameters.num_u_t,1);   
    
    //Objective function initialization
        casadi::MX J = 0;
        double current_t = 0.0;
    //Constraint function and bounds
        std::vector<casadi::MX> g;
    //initializing error vector    
        casadi::MX error = casadi::MX::sym("error",(m_model_parameters.num_x_t,1)); 
        std::vector<double> params; 
    //full trajectory V
        std::vector<casadi::MX> V_all;
        std::vector<double> v_min,v_max,v_init;

    
 
//Loop for each local setpoint
   for(int i = 0; i<15; i++ ){
    //Determine time span and shooting nodes required for setpoint
        double time_span = m_model_parameters.waypoint_list[i+1][0] - m_model_parameters.waypoint_list[i][0];
        double num_shooting_nodes_t = time_span/m_model_parameters.step_size_t.as_seconds();
    //Initializing initial and terminal state constraints and initial guess
        std::vector<double> x0_min = m_model_parameters.waypoint_list[i];
        std::vector<double> xf_min = m_model_parameters.waypoint_list[i+1];
    //Erase timestamp from setpoint
        x0_min.erase(x0_min.begin());
        xf_min.erase(xf_min.begin());

    //These likely could/ should be pointers for computation efficiency
        std::vector<double> x_init = x0_min;
        std::vector<double> x0_max = x0_min;
        std::vector<double> xf_max = xf_min;

    //Total number of variables (number of states * (time steps+1)  +  number of control inputs * time steps 
        int NV = m_model_parameters.num_x_t*(num_shooting_nodes_t+1) + m_model_parameters.num_u_t*num_shooting_nodes_t;      
    //Declare a variable vector to use in NLP
        casadi::MX V = casadi::MX::sym("V",NV);
    //NLP variable bounds and initial guesses
        
    //Offset in V --- this is like a counter variable
        int offset=0;

    //Declare vectors for the state and control at each node, insert state and control bounds
        std::vector<casadi::MX> X, U;
    //Loop through shooting nodes from initial state to desired states
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

        V_all.push_back(V);

    //Loop over shooting nodes for calculating dynamics and cost function
        for(int k=0; k<num_shooting_nodes_t; ++k){
        //Defining current state from xi+1 = f(xi,ui)     
            casadi::MX current_state = casadi::MX::sym("current_state", m_model_parameters.num_x_t);
            auto funcOut =\
            F({{"x",X[k]},{"u",U[k]}});
            current_state = funcOut["x_next"];
    // Save dynamic constraints
            g.push_back(current_state-X[k+1]);

            auto desired_state  = reshape(params_sym(casadi::Slice(0,start_Q)),m_model_parameters.num_x_t,1);
            error = current_state - desired_state;
            auto error_func_res = err_fun({{"curr_pos_fun",current_state},{"des_pos_fun",desired_state}});
            error = error_func_res["err_out"];
            auto delta_U = U[k] - ((k == 0)? u_init_in : U[k-1]); //conditional operator
        
//cost function contributions       
            if (std::count(m_cf_params.begin(),m_cf_params.end(),0)){
                J += mtimes(error.T(),mtimes(Q,error));    // Q
            }
            if (std::count(m_cf_params.begin(),m_cf_params.end(),1)){
                J += mtimes(delta_U.T(),mtimes(R,delta_U));  // R
            }
            if (std::count(m_cf_params.begin(),m_cf_params.end(),2)){
                J += mtimes(U[k].T(),mtimes(Rm,U[k])); //Rm
            }  
        }

        for(int k = 1; k<m_model_parameters.waypoint_list[i].size(); k++){
            params.push_back(m_model_parameters.waypoint_list[i+1][k]);
        }  
//Increasing traj_size based on cost function parameters
        if (std::count(m_cf_params.begin(),m_cf_params.end(),0)){
            for(int i = 0; i<Q_.size(); i++){
                params.push_back(Q_[i]);    // Q
            }
        }
        if (std::count(m_cf_params.begin(),m_cf_params.end(),1)){
            for(int i = 0; i<R_.size(); i++){
                params.push_back(R_[i]);    // Q
            } // R
        }
        if (std::count(m_cf_params.begin(),m_cf_params.end(),2)){
            for(int i = 0; i<Rm_.size(); i++){
                params.push_back(Rm_[i]);    // Q
            }  // Rm
        }

        for(int i = 0; i<u_init.size(); i++){
                params.push_back(u_init[i]);   // Q
        } 
           
    } // Waypoint list size - Number of IPOPT Iterations
//THIS IS WHERE I WANT THE WAYPOINT LOOP TO END        
        
    //Vertically concatenating dynamic constraints
        casadi::MX g_vec = casadi::MX::vertcat(g);
        casadi::MX v_all_vec = casadi::MX::vertcat(V_all);
    // NLP formulation
        casadi::MXDict nlp = {{"x", v_all_vec}, {"f", J}, {"g", g_vec}, {"p", params_sym}};
    //Default opts for IPOPT
        casadi::Dict opts;
    //Create an NLP solver and buffers
        m_solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
   
    //Adding state and control bounds and initial guesses to problem
        std::map<std::string, DM> arg, res;
        arg["p"] = params;
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
        // for(int k = 0; k<num_shooting_nodes_t;k++){
        //     data_line.clear();
        //     //states at node
        //     data_line.push_back(m_model_parameters.step_size_t.as_seconds()*k+ m_model_parameters.waypoint_list[i][0]);
        //     for (int i = 0; i < m_model_parameters.num_x_t; i++) data_line.push_back(mahi::util::RAD2DEG*V_opt[(k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i)]);//[) ;
        //     //control at node
        //     for (int i = 0; i < m_model_parameters.num_u_t; i++) data_line.push_back(V_opt[k*(m_model_parameters.num_x_t+m_model_parameters.num_u_t)+i+m_model_parameters.num_x_t]);
        //     data.push_back(data_line);   
        // }
        

        // std::vector<double> data_last_timestep= data.back();
        // std::vector<double> previous_u;
        // for (int i = data_last_timestep.size()-1; i >= (data_last_timestep.size()-m_model_parameters.num_u_t); i--) previous_u.push_back(data_last_timestep[i]);
        // reverse(previous_u.begin(),previous_u.end());
        // u_init = previous_u;
 

//Data saving         
    std::vector<std::string> header_pos_names = {"EFE ref (deg)", "FPS ref (deg)", "WFE ref (deg)", "WRU ref (deg)"};
    std::vector<std::string> header_vel_names = {"EFE ref (deg/s)", "FPS ref (deg/s)", "WFE ref (deg/s)", "WRU ref (deg/s)"};
    std::vector<std::string> header_fes_names = {"Alpha 0", "Alpha 1", "Alpha 2", "Alpha 3", "Alpha 4", "Alpha 5", "Alpha 6", "Alpha 7"};
        
    std::vector<std::string> header;
    header.push_back("Time (s)");

    for(int k = 0; k<dof.size(); k++) header.push_back(header_pos_names[dof[k]]);
    for(int k = 0; k<dof.size(); k++) header.push_back(header_vel_names[dof[k]]);
        
    for(int k = 0; k<=m_model_parameters.muscles_enabled.size(); k++){ 
        if(m_model_parameters.muscles_enabled[k]){
            header.push_back(header_fes_names[k]);
        }
    }
    std::cout<<Q_in<<std::endl;
    std::cout<<R_in<<std::endl;
    std::cout<<Rm_in<<std::endl;
    mahi::util::Timestamp ts;
    std::string save_filepath = "C:/Git/fes-exo-traj-opt/deidentified_data/" + m_model_parameters.name_t + "/Trajectories/" + dof_string + "_optimized_trajectory/" + cf_string + "/"+ std::to_string(trial_num) +"_"+ ts.yyyy_mm_dd_hh_mm_ss();
    std::cout << save_filepath << std::endl;
    mahi::util::csv_write_row(save_filepath + ".csv",header); 
    mahi::util::csv_append_rows(save_filepath + ".csv",data);     
} // Create trajectory
