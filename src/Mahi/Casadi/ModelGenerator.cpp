#include <Mahi/Casadi/ModelGenerator.hpp>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

ModelGenerator::ModelGenerator(ModelParameters model_parameters, casadi::SX x, casadi::SX x_dot, casadi::SX u) : 
    m_model_parameters(model_parameters),
    m_x(x),
    m_x_dot(x_dot),
    m_u(u)
{
    
}

ModelGenerator::~ModelGenerator(){

}

void ModelGenerator::create_model(){

    mahi::util::print("generating a {}linear model with {} shooting nodes over {} seconds with {} states, and {} control variables", (m_model_parameters.is_linear ? "" : "non"), m_model_parameters.num_shooting_nodes, m_model_parameters.timespan.as_seconds(), m_model_parameters.num_x, m_model_parameters.num_u);

    casadi::SX curr_pos_fun = SX::sym("curr_pos_fun", m_model_parameters.num_x);
    casadi::SX des_pos_fun = SX::sym("des_pos_fun", m_model_parameters.num_x);
    casadi::SX err_out = curr_pos_fun - des_pos_fun; 
    casadi::Function err_fun = casadi::Function("err_fun",{curr_pos_fun, des_pos_fun},{err_out},{"curr_pos_fun", "des_pos_fun"},{"err_out"});

    // nonlinear version
    casadi::SX x_next = m_x + m_x_dot*m_model_parameters.step_size.as_seconds();
    casadi::Function F = casadi::Function("F",{m_x,m_u},{x_next},{"x","u"},{"x_next"});

    // I made all of these symbolic so that I could pass them in as arguments, but not sure that is helpful.
    // I think these could be consolidated and cleaned up. Not sure how much benefit that would provide
    casadi::SX A_sym = casadi::SX::sym("A_sym",m_model_parameters.num_x,m_model_parameters.num_x);
    casadi::SX B_sym = casadi::SX::sym("B_sym",m_model_parameters.num_x,m_model_parameters.num_u);
    casadi::SX x_init_sym = casadi::SX::sym("x_init_sym",m_model_parameters.num_x);
    casadi::SX u_init_sym = casadi::SX::sym("u_init_sym",m_model_parameters.num_u);
    casadi::SX x_dot_init_sym = casadi::SX::sym("x_dot_init_sym",m_model_parameters.num_x);
    
    // linear functions
    casadi::SX A = jacobian(m_x_dot,m_x);
    casadi::SX B = jacobian(m_x_dot,m_u);
    casadi::SX x_dot_lin = mtimes(A_sym,(m_x-x_init_sym)) + mtimes(B_sym,(m_u-u_init_sym)) + x_dot_init_sym;
    casadi::SX x_next_lin = m_x + x_dot_lin*m_model_parameters.step_size.as_seconds();
    
    // actual functions these first three get called only once with initial conditions, then used in F_lin later
    casadi::Function get_A = casadi::Function(m_model_parameters.name + "_get_A",{m_x,m_u},{A},{"x","u"},{"A"});
    casadi::Function get_B = casadi::Function(m_model_parameters.name + "_get_B",{m_x,m_u},{B},{"x","u"},{"B"});
    casadi::Function get_x_dot_init = casadi::Function(m_model_parameters.name + "_get_x_dot_init",{m_x,m_u},{m_x_dot},{"x","u"},{"x_dot_init"});

    generate_linear_functions(get_A, get_B, get_x_dot_init);

    // this is the regular update, and only x and u are the current timestep. I think could be consolidated
    casadi::Function F_lin = casadi::Function("F_lin",{A_sym,B_sym,m_x,m_u,x_dot_init_sym,x_init_sym,u_init_sym},{x_next_lin},{"A","B","x","u","x_dot_init","x_init","u_init"},{"x_next_lin"});

    // total number of variables (number of states * (time steps+1)  +  number of control inputs * time steps )
    int NV = m_model_parameters.num_x*(m_model_parameters.num_shooting_nodes+1) + m_model_parameters.num_u*m_model_parameters.num_shooting_nodes;

    // initial guess for u and x. We are guessing that all states as well as control inputs start at 0
    std::vector<double> u_init(m_model_parameters.num_u,0.0);
    std::vector<double> x_init(m_model_parameters.num_x,0.0);

    // inital condition for state. This is essentially saying that we start at state of 0 for both position and velocity.
    std::vector<double> x0_min(m_model_parameters.num_x,0.0);
    std::vector<double> x0_max(m_model_parameters.num_x,0.0);

    // we don't set a found on final condition
    std::vector<double> xf_min(m_model_parameters.num_x,-casadi::inf);
    std::vector<double> xf_max(m_model_parameters.num_x,casadi::inf);

    // declare a variable vector to use in NLP
    casadi::MX V = casadi::MX::sym("V",NV);

    // NLP variable bounds and initial guesses
    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    // declare vectors for the state and control at each node
    std::vector<casadi::MX> X, U;
    for(int k=0; k<m_model_parameters.num_shooting_nodes; ++k){
        // Local state
        X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x))));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), m_model_parameters.x_min.begin(), m_model_parameters.x_min.end());
            v_max.insert(v_max.end(), m_model_parameters.x_max.begin(), m_model_parameters.x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += m_model_parameters.num_x;

        // Local Control
        U.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_u))));
        v_min.insert(v_min.end(), m_model_parameters.u_min.begin(), m_model_parameters.u_min.end());
        v_max.insert(v_max.end(), m_model_parameters.u_max.begin(), m_model_parameters.u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
        offset += m_model_parameters.num_u;
    }

    // State at end
    X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_model_parameters.num_x))));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += m_model_parameters.num_x;

    // Objective function
    casadi::MX J = 0;
    double current_t = 0.0;

    std::vector<casadi::MX> error_vec;
    // error_vec.resize(num_shooting_nodes)

    casadi::MX error = casadi::MX::sym("error",(m_model_parameters.num_x,1));
    casadi::MX integral = casadi::MX::sym("integral",(m_model_parameters.num_x,1));
    casadi::MX Q = casadi::MX::eye(m_model_parameters.num_x); // weight on positions and velocities
    casadi::MX P = casadi::MX::eye(m_model_parameters.num_x/2); // weight on integral term
    casadi::MX R = casadi::MX::eye(m_model_parameters.num_u); // weight on change in control inputs
    casadi::MX Rm = casadi::MX::eye(m_model_parameters.num_u);// weight on maginutde of control inputs
    
    // Constratin function and bounds
    std::vector<casadi::MX> g_vec;

    int traj_size = m_model_parameters.num_shooting_nodes*m_model_parameters.num_x;

    traj_size += m_model_parameters.num_x  // Q
               + m_model_parameters.num_u  // R
               + m_model_parameters.num_u  // Rm
               + m_model_parameters.num_x/2  // P
               + m_model_parameters.num_x/2; // integral


    if (m_model_parameters.is_linear){
        traj_size += m_model_parameters.num_x*m_model_parameters.num_x // A
                   + m_model_parameters.num_x*m_model_parameters.num_u // B
                   + m_model_parameters.num_x         // x_dot_init
                   + m_model_parameters.num_x         // x_init
                   + m_model_parameters.num_u;        // u_init
    }

    // create a vector of symbolic variables that we will import. This is of size (num_time_steps * num_states)
    // fir the nonlinear case, and for the linear case, A, B, x_dot_init, x_init, and u_init are also passed
    casadi::MX traj = casadi::MX::sym("traj", traj_size); // traj_size

    int start_Q  = (int)(m_model_parameters.num_shooting_nodes*m_model_parameters.num_x);
    int start_R  = start_Q   + m_model_parameters.num_x;
    int start_Rm = start_R   + m_model_parameters.num_u;
    int start_P  = start_Rm  + m_model_parameters.num_u;
    int start_int= start_P   + m_model_parameters.num_x/2;
    int end_int  = start_int + m_model_parameters.num_x/2; 

    casadi::MX Q_in  = reshape(traj(casadi::Slice(start_Q,start_R)),m_model_parameters.num_x,1);
    casadi::MX R_in  = reshape(traj(casadi::Slice(start_R,start_Rm)),m_model_parameters.num_u,1);
    casadi::MX Rm_in = reshape(traj(casadi::Slice(start_Rm,start_P)),m_model_parameters.num_u,1);
    casadi::MX P_in  = reshape(traj(casadi::Slice(start_P,start_int)),m_model_parameters.num_x/2,1);
    casadi::MX int_in = reshape(traj(casadi::Slice(start_int,end_int)),m_model_parameters.num_x/2,1);
    // for (size_t i = 0; i < m_model_parameters.num_x/2; i++) int_in(i) = int_in_traj(i);
    for (size_t i = 0; i < m_model_parameters.num_x; i++) Q(i,i)   = Q_in(i);
    for (size_t i = 0; i < m_model_parameters.num_x/2; i++) P(i,i) = P_in(i);
    // for (size_t i = m_model_parameters.num_x/2; i < m_model_parameters.num_x; i++) P(i,i)  = 0.0;
    for (size_t i = 0; i < m_model_parameters.num_u; i++) R(i,i)   = R_in(i);
    for (size_t i = 0; i < m_model_parameters.num_u; i++) Rm(i,i)  = Rm_in(i);

    casadi::MX lin_A;
    casadi::MX lin_B;
    casadi::MX x_dot_init;
    casadi::MX x_init_in;
    casadi::MX u_init_in;

    // if it is linear, we need to get the extra variables that were passed as parameters and rshape them
    // to the proper size to use
    if (m_model_parameters.is_linear){
        int start_A          = end_int;
        int start_B          = start_A + m_model_parameters.num_x*m_model_parameters.num_x;
        int start_x_dot_init = start_B + m_model_parameters.num_x*m_model_parameters.num_u;
        int start_x_init     = start_x_dot_init + m_model_parameters.num_x;
        int start_u_init     = start_x_init + m_model_parameters.num_x;
        int end_u_init       = start_u_init + m_model_parameters.num_u;

        lin_A      = reshape(traj(casadi::Slice(start_A,start_B)),m_model_parameters.num_x,m_model_parameters.num_x);
        lin_B      = reshape(traj(casadi::Slice(start_B,start_x_dot_init)),m_model_parameters.num_x,m_model_parameters.num_u);
        x_dot_init = reshape(traj(casadi::Slice(start_x_dot_init,start_x_init)),m_model_parameters.num_x,1);
        x_init_in  = reshape(traj(casadi::Slice(start_x_init,start_u_init)),m_model_parameters.num_x,1);
        u_init_in  = reshape(traj(casadi::Slice(start_u_init,end_u_init)),m_model_parameters.num_u,1);
    }

    // Loop over shooting nodes
    for(int k=0; k<m_model_parameters.num_shooting_nodes; ++k){
        // iterate through dynamics using either linearized or nonlinear dynamics
        // casadi::MX current_state;
        casadi::MX current_state = casadi::MX::sym("current_state", m_model_parameters.num_x); // traj_size
        if (m_model_parameters.is_linear){
            auto funcOut =\
            F_lin({{"A",lin_A},{"B",lin_B},{"x",X[k]},{"u",U[k]},{"x_dot_init",x_dot_init},{"x_init",x_init_in},{"u_init",u_init_in}}); // make these inputs instead of variables
            current_state = funcOut["x_next_lin"];
        }
        else{
            auto funcOut =\
            F({{"x",X[k]},{"u",U[k]}});
            current_state = funcOut["x_next"];
        }
        
        // Save continuity constraints
        g_vec.push_back(current_state-X[k+1]);
        auto desired_state = traj.nz(casadi::Slice(k*static_cast<int>(m_model_parameters.num_x),(k+1)*static_cast<int>(m_model_parameters.num_x)));

        // error = current_state - desired_state;
        auto error_func_res = err_fun({{"curr_pos_fun",current_state},{"des_pos_fun",desired_state}});
        error = error_func_res["err_out"];

        
        // auto integral_add = error.nz(casadi::Slice(0,m_model_parameters.num_x/2));
        if (k == 0) {
            integral = int_in + error.nz(casadi::Slice(0,static_cast<int>(m_model_parameters.num_x)/2));
        }
        else {
            if (k < 10) integral = 0.95*integral + error.nz(casadi::Slice(0,static_cast<int>(m_model_parameters.num_x)/2))*m_model_parameters.step_size.as_seconds();//pos_error
        }

        J += mtimes(error.T(),mtimes(Q,error));//error.T()*Q*error;

        // Add objective function contribution on change of input
        auto delta_U = U[k] - ((k == 0) ? u_init_in : U[k-1]);
        J += mtimes(delta_U.T(),mtimes(R,delta_U));

        // Add objective function contribution on change of input
        J += mtimes(U[k].T(),mtimes(Rm,U[k]));

        // auto result_int = mtimes(integral.T(),mtimes(P,integral));
        // auto result_err = mtimes(error.T(),mtimes(Q,error));
        // if (k == 0) std::cout << result_int << std::endl;
        // if (k == 0) std::cout << result_err << std::endl;
        // J += mtimes(integral.T(),mtimes(P,integral));
    }

    // NLP
    casadi::MX g_vec2 = casadi::MX::vertcat(g_vec);
    casadi::MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec2}, {"p", traj}};

    // leave default but can add later at runtime
    casadi::Dict opts;

    // Create an NLP solver and buffers
    m_solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
}

void ModelGenerator::generate_c_code(){
    m_c_file_filepath = m_model_parameters.name + ".c";
    std::cout << "generating c code to filepath " << m_c_file_filepath << std::endl;
    m_solver.generate_dependencies(m_c_file_filepath);
}

void ModelGenerator::generate_linear_functions(casadi::Function get_A, casadi::Function get_B, casadi::Function get_x_dot_init){
    auto C = CodeGenerator(m_model_parameters.name + "_linear_functions.c");
    C.add(get_A);
    C.add(get_B);
    C.add(get_x_dot_init);
    C.generate();
    int flag = system(("gcc -fPIC -shared -O3 "
                      + m_model_parameters.name + "_linear_functions.c"
                      + " -o " + m_model_parameters.name + "_linear_functions.so").c_str());
    casadi_assert(flag==0, "Compilation failed");
}

void ModelGenerator::compile_model(){
    m_model_parameters.dll_filepath = m_model_parameters.name + ".so";
    std::cout << "generating dll at filepath " << m_model_parameters.dll_filepath << std::endl;
    int flag = system(("gcc -fPIC -shared -O1 " + m_c_file_filepath + " -o " + m_model_parameters.dll_filepath).c_str());
    casadi_assert(flag==0, "Compilation failed");
    save_param_file();
}

void ModelGenerator::save_param_file(){
    nlohmann::json j;

    j["model"] = m_model_parameters;

    std::ofstream file1(m_model_parameters.name + ".json");
    if (file1.is_open())
        file1 << j;
    file1.close();
}
