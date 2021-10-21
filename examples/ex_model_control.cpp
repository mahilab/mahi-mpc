#include <Mahi/Casadi/ModelGenerator.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using mahi::util::PI;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("d,dll", "Generate new dll.")
        ("m,ma27", "Runs using MA27 solver. Mumps otherwise")
        ("l,linear", "Generates linearized model.");
    
    auto result = options.parse(argc, argv);

    bool linear = result.count("linear") > 0;

    // Bounds and initial guess
    casadi::Dict opts;

    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 200;
    opts["ipopt.linear_solver"] = result.count("m") ? "ma27" : "mumps";
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = 0;
    opts["ipopt.sb"] = "yes";

    std::string dll = (linear ? "linear_" : "nonlinear_") + std::string("double_pendulum.so");
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
    
    auto sim_time = mahi::util::seconds(0.2);
    double curr_sim_time = 0;
    
    std::vector<double> qA_opt,qB_opt,qA_dot_opt,qB_dot_opt;
    std::vector<double> TA_opt, TB_opt;

    double dt = time_step.as_seconds();

    double avg_solve_time = 0;
    mahi::util::Clock sim_clock;

    std::vector<double> x_next_vec = {0,0,0,0};
    std::vector<double> u_curr = {0,0};


    while (curr_sim_time < sim_time.as_seconds()){
        
        // update trajectory with new desired
        std::vector<double> traj;
        double curr_traj_time = curr_sim_time;
        for (size_t i = 0; i < ns; i++)
        {
            traj.push_back(sin(2*PI*curr_traj_time));
            traj.push_back(-sin(2*PI*curr_traj_time));
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

        static int iterations = 0;
        avg_solve_time = (avg_solve_time*iterations+(double)solve_time)/(iterations+1);
        iterations++;
        
        std::vector<double> V_opt(res.at("x"));
        arg["x0"] = V_opt;
        
        std::vector<double> x_curr = std::vector<double>(V_opt.begin(),V_opt.begin()+nx); // first 4 indices
        u_curr = std::vector<double>(V_opt.begin()+(nx),V_opt.begin()+ (nx + nu)); // 5-6 indices

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
        
        qA_opt.push_back(x_curr[0]);
        qB_opt.push_back(x_curr[1]);
        qA_dot_opt.push_back(x_curr[2]);
        qB_dot_opt.push_back(x_curr[3]);

        TA_opt.push_back(u_curr[0]);
        TB_opt.push_back(u_curr[1]);

        curr_sim_time += time_step.as_seconds();
    }
    mahi::util::print("sim time: {:.2f} s, avg solve time: {:.2f} ms",sim_clock.get_elapsed_time().as_seconds(),avg_solve_time);

    mahi::util::disable_realtime();

    // Extract the optimal state trajectory
    std::ofstream file;
    std::string filename = (linear ? "" : "non") + std::string("linear_double_pendulum_mpc_results.m");
    file.open(filename.c_str());
    file << "% Results from " __FILE__ << std::endl;
    file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
    file << std::endl;

    // Save results
    file << "t = linspace(0," << sim_time.as_seconds() << "," << (sim_time/time_step) << ");" << std::endl;
    file << "qA = " << qA_opt << ";" << std::endl;
    file << "qB = " << qB_opt << ";" << std::endl;
    file << "qA_dot = " << qA_dot_opt << ";" << std::endl;
    file << "qB_dot = " << qB_dot_opt << ";" << std::endl;
    file << "TA = " << TA_opt << ";" << std::endl;
    file << "TB = " << TB_opt << ";" << std::endl;

    file << "figure;" << std::endl;
    file << "hold on;" << std::endl;
    file << "plot(t,qA);" << std::endl;
    file << "plot(t,qB);" << std::endl;
    file << "xlabel('Time (s)');" << std::endl;
    file << "ylabel('Position (rad)');" << std::endl;
    file << "legend('qA','qB');" << std::endl; 
    std::cout << "finished" << std::endl;

    return 0;
}
