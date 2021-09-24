#include <casadi/casadi.hpp>

int main(){
    double L = 1.0;
    double m = 1.0;
    double g = 9.81;

    casadi::SX qA = casadi::SX::sym("qA");
    casadi::SX qB = casadi::SX::sym("qB");
    casadi::SX qA_dot = casadi::SX::sym("qA_dot");
    casadi::SX qB_dot = casadi::SX::sym("qB_dot");
    casadi::SX TA = casadi::SX::sym("TA");
    casadi::SX TB = casadi::SX::sym("TB");
    // state vector
    casadi::SX x = casadi::SX::vertcat({qA,qB,qA_dot,qB_dot});
    // control vector
    casadi::SX u = casadi::SX::vertcat({TA,TB});
    // size of state
    int nx = x.size1();
    // size of control
    int nu = u.size1();

    // Bounds for control
    std::vector<double> u_min = {-1.0*(casadi::inf), -1.0*(casadi::inf)};
    std::vector<double> u_max = {casadi::inf, casadi::inf};
    // Initial Guess for control
    std::vector<double> u_init = {0.0, 0.0};

    // Bounds on initial state
    std::vector<double> x0_min = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> x0_max = {0.0, 0.0, 0.0, 0.0}; // min and max are the same here
    // Bounds on state
    std::vector<double> x_min = {-1.0*(casadi::inf), -1.0*(casadi::inf), -1.0*(casadi::inf), -1.0*(casadi::inf)};
    std::vector<double> x_max = {casadi::inf, casadi::inf, casadi::inf, casadi::inf};
    // Bounds on final state - These could be interesting to adjust
    std::vector<double> xf_min = {-1.0*(casadi::inf), -1.0*(casadi::inf), -1.0*(casadi::inf), -1.0*(casadi::inf)};
    std::vector<double> xf_max = {casadi::inf, casadi::inf, casadi::inf, casadi::inf};
    // Initial guess for state
    std::vector<double> x_init = {0.0, 0.0, 0.0, 0.0};

    // Final Time
    double tf = 10.0;

    // Number of shooting nodes
    int ns = 50;

    // ODE right hand side
    casadi::SX qA_ddot = -(TA - TB - TB*cos(qB) + L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(L*L*m*(pow(cos(qB),2) - 2));
    casadi::SX qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + 2*L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + L*L*m*qB_dot*qB_dot*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));

    casadi::SX ode = casadi::SX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});

    // Cost function? In the example they call this quadrature
    casadi::SX quad = qA*qA + qB*qB + TA*TA + TB*TB;

    // create the DAE
    casadi::SXDict dae = {{"x",x},{"p",u},{"ode",ode},{"quad",quad}};

    // Create an integrator (this example used cvodes, but that could change)
    casadi::Function F = casadi::integrator("integrator","cvodes",dae,{{"t0",0},{"tf",tf/ns}});

    // Total number of NLP variables
    int NV = nx*(ns+1) + nu*ns;

    // declare a variable vector to use in NLP
    casadi::MX V = casadi::MX::sym("V",NV);

    // NLP variable bounds and initial guesses
    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset = 0;

    // declare vectors for the state and control at each node
    std::vector<casadi::MX> X,U;
    for(int k=0; k<ns; ++k){
        // Local state
        X.push_back(V.nz(casadi::Slice(offset,offset+nx)));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.end(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += nx;

        // Local Control
        U.push_back(V.nz(casadi::Slice(offset,offset+nu)));
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
        offset += nu;
    }

    // State at end
    X.push_back(V.nz(casadi::Slice(offset,offset+nx)));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += nx;

    // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
    casadi_assert(offset==NV, "");

    // Objective function
    casadi::MX J = 0;

    // Constratin function and bounds
    std::vector<casadi::MX> g_vec;
    // Loop over shooting nodes
    for(int k=0; k<ns; ++k){
        // Create an evaluation node
        casadi::MXDict I_out = F(casadi::MXDict{{"x0", X[k]},{"p", U[k]}});

        // Save continuity constraints
        g_vec.push_back(I_out.at("xf") - X[k+1]);

        // Add objective function contribution
        J += I_out.at("qf");
    }

    // NLP
    casadi::MX g_vec2 = casadi::MX::vertcat(g_vec);
    casadi::MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec2}};

    // Set Options
    casadi::Dict opts;
    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 200;

    // Create an NLP solver and buffers
    casadi::Function solver = casadi::nlpsol("nlpsol", "ipopt", nlp, opts);
    std::map<std::string, casadi::DM> arg, res;

    // Bounds and initial guess
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;

    std::cout << "display" << std::endl;
    
    // Solve the problem
    res = solver(arg);

    std::cout << "finished" << std::endl;
    return 0;
}