#include <Mahi/Casadi/M.hpp>
#include <Mahi/Casadi/G.hpp>
#include <Mahi/Casadi/V.hpp>
#include <casadi/casadi.hpp>
#include <Mahi/Casadi/ModelGenerator.hpp>

using namespace casadi;
using mahi::util::PI;

int main(int argc, char *argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("d,dll", "Generate new dll.")
        ("m,ma27", "Runs using MA27 solver. Mumps otherwise")
        ("l,linear", "Generates linearized model.");

    auto result = options.parse(argc, argv);

    bool linear = result.count("linear") > 0;

    const int num_x = 8;
    const int num_u = 4;

    SX q0 = SX::sym("q0");
    SX q1 = SX::sym("q1");
    SX q2 = SX::sym("q2");
    SX q3 = SX::sym("q3");
    SX q0_dot = SX::sym("q0_dot");
    SX q1_dot = SX::sym("q1_dot");
    SX q2_dot = SX::sym("q2_dot");
    SX q3_dot = SX::sym("q3_dot");
    
    SX T0 = SX::sym("T0");
    SX T1 = SX::sym("T1");
    SX T2 = SX::sym("T2");
    SX T3 = SX::sym("T3");

    // state vector
    SX x = SX::vertcat({q0,q1,q2,q3,q0_dot,q1_dot,q2_dot,q3_dot});
    SX u = SX::vertcat({T0,T1,T2,T3});

    SX q_dot = SX::vertcat({q0_dot,q1_dot,q2_dot,q3_dot});

    std::vector<double>  B_coef = {0.1215, 0.0252, 0.0019, 0.0029};
    std::vector<double> Fk_coef = {   0.5, 0.1891, 0.0541, 0.1339};

    SX B = SX::vertcat({B_coef[0]*q0_dot*1.0, B_coef[1]*q1_dot*1.0, B_coef[2]*q2_dot*1.0, B_coef[3]*q3_dot*1.0});
    SX Fk = SX::vertcat({Fk_coef[0]*tanh(q0_dot*10.0), Fk_coef[1]*tanh(q1_dot*10.0), Fk_coef[2]*tanh(q2_dot*10.0), Fk_coef[3]*tanh(q3_dot*10.0)});

    auto V = get_V(x);
    auto G = get_G(x);
    auto B_eom = u - mtimes(V,q_dot) - G - B - Fk;
    auto A_eom = get_M(x);
    SX q_d_dot = solve(A_eom,B_eom);
    SX x_dot = vertcat(q_dot,q_d_dot);

    auto A = jacobian(x_dot,x);
    auto B = jacobian(x_dot,u);
    casadi::Function get_A = casadi::Function("get_A",{x,u},{A},{"x","u"},{"A"});
    casadi::Function get_B = casadi::Function("get_B",{x,u},{B},{"x","u"},{"B"});
    casadi::Function get_x_dot_init = casadi::Function("F_get_x_dot_init",{x,u},{x_dot},{"x","u"},{"x_dot_init"});

    // Bounds on state
    std::vector<double> x_min(num_x,-inf);
    std::vector<double> x_max(num_x, inf);

    // Bounds for control
    std::vector<double> u_min(num_u,-inf);
    std::vector<double> u_max(num_u, inf);

    mahi::util::Time final_time = mahi::util::milliseconds(100);
    mahi::util::Time time_step  = mahi::util::milliseconds(5);

    ModelGenerator my_generator("moe_model", x, x_dot, u, final_time, time_step, u_min, u_max, x_min, x_max, linear);

    my_generator.create_model();
    my_generator.generate_c_code("moe_model.c");
    if (result.count("dll")) my_generator.compile_model("moe_model.so");

    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    int ns = (final_time/time_step);

    int nx = x.size1();
    int nu = u.size1();
    mahi::util::print("{} shooting nodes over {} seconds with {} states, and {} control variables", ns, final_time, nx, nu);

    // Bounds on initial state
    std::vector<double> x0_min(num_x,0);
    std::vector<double> x0_max(num_x,0); // min and max are the same here
    std::vector<double> x_init(num_x,0);
    std::vector<double> u_init(num_u,0);
    std::vector<double> xf_min(num_x,-inf);
    std::vector<double> xf_max(num_x, inf);


    //declare vectors for the state and control at each node
    for(int k=0; k<ns; ++k){
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());

        // Local Control
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    }

    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());

    // Bounds and initial guess
    casadi::Dict opts;

    std::string linear_solver = result.count("m") ? "ma27" : "mumps";

    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 200;
    opts["ipopt.linear_solver"] = linear_solver;
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = 0;
    opts["ipopt.sb"] = "yes";

    mahi::util::print("Using solver: {}",linear_solver);

    std::string dll = (linear ? "linear_" : "nonlinear_") + std::string("moe_model.so");

    Function solver = nlpsol("nlpsol", "ipopt", dll, opts);

    std::map<std::string, casadi::DM> arg, res;
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;
    
    casadi::SX x_next_euler = x + x_dot*time_step.as_seconds();
    casadi::Function F_euler = casadi::Function("F",{x,u},{x_next_euler},{"x","u"},{"x_next_euler"});

    casadi::Function x_dot_fun("x_dot_ode",{x,u},{x_dot},{"x","u"},{"x_dot"});
    
    auto sim_time = mahi::util::seconds(1);
    double curr_sim_time = 0;
    
    std::vector<double> q0_opt,q1_opt,q2_opt,q3_opt,q0_dot_opt,q1_dot_opt,q2_dot_opt,q3_dot_opt;
    std::vector<double> T0_opt, T1_opt, T2_opt, T3_opt;

    double dt = time_step.as_seconds();

    double avg_solve_time = 0;
    mahi::util::Clock sim_clock;

    std::vector<double> x_next_vec(8,0);
    std::vector<double> u_curr(4,0);

    while (curr_sim_time < sim_time.as_seconds()){
        // std::cout << "running iteration " << iter++ << std::endl;
        // update trajectory with new desired
        std::vector<double> traj;
        double curr_traj_time = curr_sim_time;
        for (size_t i = 0; i < ns; i++)
        {
            traj.push_back(sin(2*PI*curr_traj_time));
            traj.push_back(-sin(2*PI*curr_traj_time));
            traj.push_back(sin(2*PI*curr_traj_time));
            traj.push_back(-sin(2*PI*curr_traj_time));
            traj.push_back(2*PI*cos(2*PI*curr_traj_time));
            traj.push_back(-2*PI*cos(2*PI*curr_traj_time));
            traj.push_back(2*PI*cos(2*PI*curr_traj_time));
            traj.push_back(-2*PI*cos(2*PI*curr_traj_time));
            curr_traj_time += time_step.as_seconds();
        }

        if (linear){
            // get and add flattened A, B, x_dot_init, x_init 
            std::vector<double> A_res(get_A(SXDict({{"x",x_next_vec},{"u",u_curr}})).at("A"));
            std::vector<double> B_res(get_B(SXDict({{"x",x_next_vec},{"u",u_curr}})).at("B"));
            std::vector<double> x_dot_init_res(get_x_dot_init(SXDict({{"x",x_next_vec},{"u",u_curr}})).at("x_dot_init"));

            for (const auto &a_ : A_res) traj.push_back(a_);
            for (const auto &b_ : B_res) traj.push_back(b_);
            for (const auto &x_dot_init_ : x_dot_init_res) traj.push_back(x_dot_init_);
            for (const auto &x_next_ : x_next_vec) traj.push_back(x_next_);
            for (const auto &u_ : u_curr) traj.push_back(u_);
        }

        arg["p"]  = traj;
        
        mahi::util::Clock my_clock;
        
        res = solver(arg);
        auto solve_time = my_clock.get_elapsed_time().as_milliseconds();
        // std::cout << solve_time << std::endl;

        static int iterations = 0;
        avg_solve_time = (avg_solve_time*iterations+(double)solve_time)/(iterations+1);
        iterations++;
        
        std::vector<double> V_opt(res.at("x"));
        arg["x0"] = V_opt;
        
        std::vector<double> x_curr(V_opt.begin(),V_opt.begin()+nx);
        u_curr = std::vector<double>(V_opt.begin()+(nx),V_opt.begin()+ (nx+nu));
        
        // begin euler
        // auto x_next = F_euler(casadi::SXDict({{"x",x_curr},{"u",u_curr}}));
        // std::vector<double> x_next_vec(x_next.at("x_next_euler"));
        // std::copy(x_next_vec.begin(),x_next_vec.end(),v_min.begin());
        // std::copy(x_next_vec.begin(),x_next_vec.end(),v_max.begin());
        // end euler

        // begin rk4
        auto k1 = x_dot_fun(casadi::SXDict({{"x",x_curr},{"u",u_curr}}));
        auto k2 = x_dot_fun(casadi::SXDict({{"x",x_curr+dt/2*k1.at("x_dot")},{"u",u_curr}}));
        auto k3 = x_dot_fun(casadi::SXDict({{"x",x_curr+dt/2*k2.at("x_dot")},{"u",u_curr}}));
        auto k4 = x_dot_fun(casadi::SXDict({{"x",x_curr+dt*k3.at("x_dot")},{"u",u_curr}}));
        x_next_vec = std::vector<double>(x_curr + dt/6*(k1.at("x_dot")+2*k2.at("x_dot")+2*k3.at("x_dot")+k4.at("x_dot")));
        // end rk4

        std::copy(x_next_vec.begin(),x_next_vec.end(),v_min.begin());
        std::copy(x_next_vec.begin(),x_next_vec.end(),v_max.begin());

        arg["lbx"] = v_min;
        arg["ubx"] = v_max;
        
        q0_opt.push_back(x_curr[0]);
        q1_opt.push_back(x_curr[1]);
        q2_opt.push_back(x_curr[2]);
        q3_opt.push_back(x_curr[3]);
        q0_dot_opt.push_back(x_curr[4]);
        q1_dot_opt.push_back(x_curr[5]);
        q2_dot_opt.push_back(x_curr[6]);
        q3_dot_opt.push_back(x_curr[7]);

        T0_opt.push_back(u_curr[0]);
        T1_opt.push_back(u_curr[1]);
        T2_opt.push_back(u_curr[2]);
        T3_opt.push_back(u_curr[3]);

        curr_sim_time += time_step.as_seconds();
    }
    // std::cout << "sim time: " << sim_clock.get_elapsed_time() << std::endl;
    mahi::util::print("sim time: {:.2f} s, avg solve time: {:.2f} ms",sim_clock.get_elapsed_time().as_seconds(),avg_solve_time);
    // vector<double> V_opt(res.at("x"));

    mahi::util::disable_realtime();
    
    // arg["x0"] = V_opt;

    // Extract the optimal state trajectory
    std::ofstream file;
    std::string filename = (linear ? "" : "non") + std::string("linear_moe_mpc_results.m");
    file.open(filename.c_str());
    file << "% Results from " __FILE__ << std::endl;
    file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
    file << std::endl;

    // Save results
    file << "t = linspace(0," << sim_time.as_seconds() << "," << (sim_time/time_step) << ");" << std::endl;
    file << "q0 = " << q0_opt << ";" << std::endl;
    file << "q1 = " << q1_opt << ";" << std::endl;
    file << "q2 = " << q2_opt << ";" << std::endl;
    file << "q3 = " << q3_opt << ";" << std::endl;
    file << "q0_dot = " << q0_dot_opt << ";" << std::endl;
    file << "q1_dot = " << q1_dot_opt << ";" << std::endl;
    file << "q2_dot = " << q2_dot_opt << ";" << std::endl;
    file << "q3_dot = " << q3_dot_opt << ";" << std::endl;
    file << "T0 = " << T0_opt << ";" << std::endl;
    file << "T1 = " << T1_opt << ";" << std::endl;
    file << "T2 = " << T2_opt << ";" << std::endl;
    file << "T3 = " << T3_opt << ";" << std::endl;

    file << "figure;" << std::endl;
    file << "hold on;" << std::endl;
    file << "plot(t,q0);" << std::endl;
    file << "plot(t,q1);" << std::endl;
    file << "plot(t,q2);" << std::endl;
    file << "plot(t,q3);" << std::endl;
    file << "xlabel('Time (s)');" << std::endl;
    file << "ylabel('Position (rad)');" << std::endl;
    file << "legend('q0','q1','q2','q3');" << std::endl; 
    std::cout << "finished" << std::endl;

    return 0;
}
