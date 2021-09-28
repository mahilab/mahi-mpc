#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using namespace std;

int main(){
    double L = 1.0;
    double m = 1.0;
    double g = 9.81;

    SX qA = SX::sym("qA");
    SX qB = SX::sym("qB");
    SX qA_dot = SX::sym("qA_dot");
    SX qB_dot = SX::sym("qB_dot");
    SX TA = SX::sym("TA");
    SX TB = SX::sym("TB");
    // state vector
    SX x = SX::vertcat({qA,qB,qA_dot,qB_dot});
    // control vector
    SX u = SX::vertcat({TA,TB});
    // size of state
    int nx = x.size1();
    // size of control
    int nu = u.size1();

    // Bounds for control
    vector<double> u_min = {-inf, -inf};
    vector<double> u_max = {inf, inf};
    // Initial Guess for control
    vector<double> u_init = {0.0, 0.0};

    // // Bounds on initial state
    vector<double> x0_min = {0, 0, 0, 0};
    vector<double> x0_max = {0, 0, 0, 0}; // min and max are the same here
    // // Bounds on state
    vector<double> x_min = {-inf, -inf, -inf, -inf};
    vector<double> x_max = {inf, inf, inf, inf};
    // // Bounds on final state - These could be interesting to adjust
    vector<double> xf_min = {-inf, -inf, -inf, -inf};
    vector<double> xf_max = {inf, inf, inf, inf};
    // // Initial guess for state
    vector<double> x_init = {0, 0, 0, 0};

    // // Final Time
    double tf = 1.0;

    // // Number of shooting nodes
    int ns = 20;

    // // ODE right hand side
    // SX qA_ddot = -(TA - TB - TB*cos(qB) + L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(L*L*m*(pow(cos(qB),2) - 2));
    // SX qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + 2*L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + L*L*m*qB_dot*qB_dot*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));

    // SX ode = SX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});

    // // Cost function? In the example they call this quadrature
    // SX quad = qA*qA + qB*qB + qA_dot*qA_dot + qB_dot*qB_dot; // + TA*TA + TB*TB;
    
    // // create the DAE
    // SXDict dae = {{"x", x}, {"p", u}, {"ode", ode}, {"quad", quad}};
    
    double dt = tf/ns;
    // SX x_next = x + ode*dt;
    // Function F = Function("F",{x,u},{x_next},{"x","u"},{"x_next"});
    // // Create an integrator (this example used cvodes, but that could change)
    // //Function F = integrator("integrator", "cvodes", dae, {{"t0", 0}, {"tf", tf/ns}});
    // //F.disp(cout);
    // // Total number of NLP variables
    int NV = nx*(ns+1) + nu*ns;

    // // declare a variable vector to use in NLP
    MX V = MX::sym("V",NV);

    // // NLP variable bounds and initial guesses
    vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    // declare vectors for the state and control at each node
    vector<MX> X, U;
    for(int k=0; k<ns; ++k){
        // Local state
        X.push_back( V.nz(Slice(offset,offset+nx)));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());
        offset += nx;

        // Local Control
        U.push_back(V.nz(Slice(offset,offset+nu)));
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
        offset += nu;
    }

    // State at end
    X.push_back(V.nz(Slice(offset,offset+nx)));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += nx;

    // // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
    // //casadi_assert(offset==NV, "");

    // // Objective function
    // MX J = 0;
    // double desired_pos = 0.0;
    // double desired_vel = 0.0;
    // double current_t = 0.0;
    
    // std::vector<double> desired_state(4,0);
    // casadi::MX error = casadi::MX::sym("error",(4,1));
    // casadi::MX Q = casadi::MX::eye(4);
    // // Constratin function and bounds
    // vector<MX> g_vec;
    // // vector<MX> current_state;
    // // Loop over shooting nodes
    // for(int k=0; k<ns; ++k){
    //     // Create an evaluation node
    //     //MXDict I_out = F(MXDict{{"x0", X[k]},{"p", U[k]}});
    //     auto funcOut =\
    //     F({{"x",X[k]},{"u",U[k]}});
    //     auto current_state = funcOut["x_next"];
    //     // Save continuity constraints
    //     //g_vec.push_back(I_out.at("xf") - X[k+1]);
    //     g_vec.push_back(current_state-X[k+1]);
    //     current_t = k*dt;
    //     desired_pos = mahi::util::PI/2*sin(current_t);
    //     desired_vel = mahi::util::PI/2*cos(current_t);
    //     desired_state = {desired_pos,desired_pos,desired_vel,desired_vel};
    //     // Add objective function contribution
    //     //J += I_out.at("qf");
    //     error = current_state - desired_state;
    //     //cout << "print here" << endl;
    //     J += mtimes(error.T(),mtimes(Q,error));//error.T()*Q*error;
    // }

    // // NLP
    // MX g_vec2 = MX::vertcat(g_vec);
    // MXDict nlp = {{"x", V}, {"f", J}, {"g", g_vec2}};

    // // Set Options
    Dict opts;
    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 200;
    // opts["ipopt.print_level"] = 0;
    opts["ipopt.linear_solver"] = "ma57";
    // opts["print_time"] = 0;
    // opts["ipopt.sb"] = "yes";

    // // Create an NLP solver and buffers
    Function solver = nlpsol("nlpsol", "ipopt", "nlp_multishoot.so", opts);
    // Function solver = nlpsol("nlpsol", "ipopt", "multishoot_library.dll", opts);

    // solver.generate_dependencies("multishoot_c.c");

    map<string, DM> arg, res;

    // Bounds and initial guess
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;

    
    mahi::util::Clock my_clock;
    // Solve the problem
    res = solver(arg);

    std::cout << my_clock.get_elapsed_time();

    

    // The optimal solution
    vector<double> V_opt(res.at("x"));
    
    arg["x0"] = V_opt;

    my_clock.restart();

    res = solver(arg);

    std::cout << my_clock.get_elapsed_time() << endl;


    // Extract the optimal state trajectory
    vector<double> qA_opt(ns+1),qB_opt(ns+1),qA_dot_opt(ns+1),qB_dot_opt(ns+1);

    for(int i=0; i<=ns; ++i){
        qA_opt[i] = V_opt.at(i*(nx+2));
        qB_opt[i] = V_opt.at(1+i*(nx+2));
        qA_dot_opt[i] = V_opt.at(2+i*(nx+2));
        qB_dot_opt[i] = V_opt.at(3+i*(nx+2));
    }

    // Get the optimal control
    vector<double> TA_opt(ns), TB_opt(ns);
    for(int i=0; i<ns; ++i){
        TA_opt[i] = V_opt.at(nx + i*(nx+2));
        TB_opt[i] = V_opt.at(nx + 1 + i*(nx+2));
    }

    ofstream file;
    string filename = "my_multishoot_results.m";
    file.open(filename.c_str());
    file << "% Results from " __FILE__ << endl;
    file << "% Generated " __DATE__ " at " __TIME__ << endl;
    file << endl;

    // Save results
    file << "t = linspace(0,10," << ns << "+1);" << endl;
    file << "qA = " << qA_opt << ";" << endl;
    file << "qB = " << qB_opt << ";" << endl;
    file << "qA_dot = " << qA_dot_opt << ";" << endl;
    file << "qB_dot = " << qB_dot_opt << ";" << endl;
    file << "TA = " << TA_opt << ";" << endl;
    file << "TB = " << TB_opt << ";" << endl;

    file << "figure;" << endl;
    file << "hold on;" << endl;
    file << "plot(t,qA);" << endl;
    file << "plot(t,qB);" << endl;
    file << "xlabel('Time (s)');" << endl;
    file << "ylabel('Position (rad)');" << endl;
    file << "legend('qA','qB');" << endl; 
    cout << "finished" << endl;
    return 0;
}