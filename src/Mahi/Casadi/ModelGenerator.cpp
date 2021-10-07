#include <Mahi/Casadi/ModelGenerator.hpp>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

ModelGenerator::ModelGenerator(std::string name, casadi::SX x, casadi::SX x_dot, casadi::SX u, mahi::util::Time final_time, mahi::util::Time step_size, std::vector<double> u_min, std::vector<double> u_max, std::vector<double> x_min, std::vector<double> x_max, bool linear) : 
    m_name(name),
    m_timespan(final_time),
    m_step_size(step_size),
    m_x(x),
    m_x_dot(x_dot),
    m_u(u),
    m_num_x(m_x.size1()),
    m_num_u(m_u.size1()),
    m_num_shooting_nodes(m_timespan.as_milliseconds()/m_step_size.as_milliseconds()),
    m_x_min(x_min),
    m_u_min(u_min),
    m_x_max(x_max),
    m_u_max(u_max),
    m_linear(linear)
{

}

ModelGenerator::~ModelGenerator(){

}

void ModelGenerator::create_model(){

    mahi::util::print("generating a {}linear model with {} shooting nodes over {} seconds with {} states, and {} control variables", (m_linear ? "" : "non"), m_num_shooting_nodes, m_timespan.as_seconds(), m_num_x, m_num_u);

    // nonlinear version
    casadi::SX x_next = m_x + m_x_dot*m_step_size.as_seconds();
    casadi::Function F = casadi::Function("F",{m_x,m_u},{x_next},{"x","u"},{"x_next"});

    // linearization
    auto A = jacobian(m_x_dot,m_x);
    auto B = jacobian(m_x_dot,m_u);
    casadi::Function get_A = casadi::Function("get_A",{m_x,m_u},{A},{"x","u"},{"A"});
    casadi::Function get_B = casadi::Function("get_B",{m_x,m_u},{B},{"x","u"},{"B"});
    casadi::Function get_x_dot_init = casadi::Function("F_get_x_dot_init",{m_x,m_u},{m_x_dot},{"x","u"},{"x_dot_init"});
    // linear functions
    casadi::SX A_sym = casadi::SX::sym("A_sym",m_num_x,m_num_x);
    casadi::SX B_sym = casadi::SX::sym("B_sym",m_num_x,m_num_u);
    casadi::SX x_init_sym = casadi::SX::sym("x_init_sym",m_num_x);
    casadi::SX u_init_sym = casadi::SX::sym("u_init_sym",m_num_u);
    casadi::SX x_dot_init_sym = casadi::SX::sym("x_dot_init_sym",m_num_x);
    auto x_dot_lin = mtimes(A_sym,(m_x-x_init_sym)) + mtimes(B_sym,(m_u-u_init_sym)) + x_dot_init_sym;
    
    casadi::SX x_next_lin = m_x + x_dot_lin*m_step_size.as_seconds();

    casadi::Function F_lin = casadi::Function("F_lin",{A_sym,B_sym,m_x,m_u,x_dot_init_sym,x_init_sym,u_init_sym},{x_next_lin},{"A","B","x","u","x_dot_init","x_init","u_init"},{"x_next_lin"});


    int NV = m_num_x*(m_num_shooting_nodes+1) + m_num_u*m_num_shooting_nodes;

    std::vector<double> u_init(m_num_u,0.0);
    std::vector<double> x_init(m_num_x,0.0);

    std::vector<double> x0_min(m_num_x,0.0);
    std::vector<double> x0_max(m_num_x,0.0);

    std::vector<double> xf_min(m_num_x,-casadi::inf);
    std::vector<double> xf_max(m_num_x,casadi::inf);

    // declare a variable vector to use in NLP
    casadi::MX V = casadi::MX::sym("V",NV);

    // NLP variable bounds and initial guesses
    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    // declare vectors for the state and control at each node
    std::vector<casadi::MX> X, U;
    for(int k=0; k<m_num_shooting_nodes; ++k){
        // Local state
        X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_num_x))));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), m_x_min.begin(), m_x_min.end());
            v_max.insert(v_max.end(), m_x_max.begin(), m_x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += m_num_x;

        // Local Control
        U.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_num_u))));
        v_min.insert(v_min.end(), m_u_min.begin(), m_u_min.end());
        v_max.insert(v_max.end(), m_u_max.begin(), m_u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
        offset += m_num_u;
    }

    // State at end
    X.push_back(V.nz(casadi::Slice(offset,offset+static_cast<int>(m_num_x))));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += m_num_x;

    // Objective function
    casadi::MX J = 0;
    double desired_pos = 0.0;
    double desired_vel = 0.0;
    double current_t = 0.0;
    
    // casadi::MX desired_state(m_num_x,0);
    casadi::MX error = casadi::MX::sym("error",(m_num_x,1));
    casadi::MX Q = casadi::MX::eye(m_num_x);

    // Constratin function and bounds
    std::vector<casadi::MX> g_vec;

    casadi::MX traj = casadi::MX::sym("traj",m_num_shooting_nodes*m_num_x);

    // get A and B based on the current values of 
    auto lin_A = get_A({{"x",X[0]},{"u",U[0]}}).at("A");
    auto lin_B = get_B({{"x",X[0]},{"u",U[0]}}).at("B");
    auto x_dot_init = get_x_dot_init({{"x",X[0]},{"u",U[0]}}).at("x_dot_init");
    // auto lin_A = lin_A_res["A"];
    // auto lin_B = lin_B_res["B"];
    // auto x_dot_init = x_dot_init_res["x_dot_init"];

    // Loop over shooting nodes
    for(int k=0; k<m_num_shooting_nodes; ++k){
        // iterate through dynamics using either linearized or nonlinear dynamics
        casadi::MX current_state;
        if (m_linear){
            auto funcOut =\
            F_lin({{"A",lin_A},{"B",lin_B},{"x",X[k]},{"u",U[k]},{"x_dot_init",x_dot_init},{"x_init",X[0]},{"u_init",U[0]}});
            current_state = funcOut["x_next_lin"];
        }
        else{
            auto funcOut =\
            F({{"x",X[k]},{"u",U[k]}});
            current_state = funcOut["x_next"];
        }
        
        // Save continuity constraints
        //g_vec.push_back(I_out.at("xf") - X[k+1]);
        g_vec.push_back(current_state-X[k+1]);
        auto desired_state = traj.nz(casadi::Slice(k*static_cast<int>(m_num_x),(k+1)*static_cast<int>(m_num_x)));
        // Add objective function contribution
        //J += I_out.at("qf");
        error = current_state - desired_state;
        J += mtimes(error.T(),mtimes(Q,error));//error.T()*Q*error;
    }

    // NLP
    casadi::MX g_vec2 = casadi::MX::vertcat(g_vec);
    casadi::MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec2}, {"p", traj}};

    // leave default but can add later at runtime
    casadi::Dict opts;

    // Create an NLP solver and buffers
    m_solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
}

void ModelGenerator::generate_c_code(std::string filepath){
    m_c_file_filepath = (m_linear ? "linear_" : "nonlinear_") + filepath;
    std::cout << "generating c code to filepath " << m_c_file_filepath << std::endl;
    m_solver.generate_dependencies(m_c_file_filepath);
}

void ModelGenerator::compile_model(std::string filepath){
    m_dll_filepath = (m_linear ? "linear_" : "nonlinear_") + filepath;
    std::cout << "generating dll at filepath " << m_dll_filepath << std::endl;
    int flag = system(("gcc -fPIC -shared -O3 " + m_c_file_filepath + " -o " + m_dll_filepath).c_str());
    casadi_assert(flag==0, "Compilation failed");
}
