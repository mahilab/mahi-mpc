#include <Mahi/Casadi/ModelGenerator.hpp>

using namespace casadi;
using mahi::util::PI;

int main(int argc, char const *argv[])
{
    double L = 1.0;
    double m = 1.0;
    double g = 9.81;

    SX qA = SX::sym("qA");
    SX qB = SX::sym("qB");
    SX qA_dot = SX::sym("qA_dot");
    SX qB_dot = SX::sym("qB_dot");
    SX TA = SX::sym("TA");
    SX TB = SX::sym("TB");

    // ODE right hand side
    SX qA_ddot = -(TA - TB - TB*cos(qB) + L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));
    SX qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + 2*L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + L*L*m*qB_dot*qB_dot*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));

    // state vector
    SX x = SX::vertcat({qA,qB,qA_dot,qB_dot});

    SX x_dot = SX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});
    
    // control vector
    SX u = SX::vertcat({TA,TB});

    // Bounds on state
    std::vector<double> x_min = {-inf, -inf, -inf, -inf};
    std::vector<double> x_max = {inf, inf, inf, inf};

    // Bounds for control
    std::vector<double> u_min = {-inf, -inf};
    std::vector<double> u_max = {inf, inf};

    mahi::util::Time final_time = mahi::util::seconds(0.5);
    mahi::util::Time time_step  = mahi::util::milliseconds(10);

    ModelGenerator my_generator("double_pendulum", x, x_dot, u, final_time, time_step, u_min, u_max, x_min, x_max);

    my_generator.create_model();
    my_generator.generate_c_code("double_pendulum.c");
    // my_generator.compile_model("double_pendulum.so");

    // std::vector<double> traj;
    // for (size_t i = 0; i < ns; i++)
    // {
    //     traj.push_back(sin(2*PI*time));
    //     traj.push_back(-sin(2*PI*time));
    //     traj.push_back(2*PI*cos(2*PI*time));
    //     traj.push_back(-2*PI*cos(2*PI*time));
    //     time += time_step.as_seconds();
    // }
        // NLP variable bounds and initial guesses
    std::vector<double> v_min,v_max,v_init;

    // Offset in V --- this is like a counter variable
    int offset=0;

    int ns = (final_time/time_step);

    int nx = x.size1();
    int nu = u.size1();
    mahi::util::print("{} shooting nodes over {} seconds with {} states, and {} control variables", ns, final_time, nx, nu);

    // Bounds on initial state
    std::vector<double> x0_min = {0, 0, 0, 0};
    std::vector<double> x0_max = {0, 0, 0, 0}; // min and max are the same here
    std::vector<double> x_init = {0, 0, 0, 0};
    std::vector<double> u_init = {0.0, 0.0};
    std::vector<double> xf_min = {-inf, -inf, -inf, -inf};
    std::vector<double> xf_max = {inf, inf, inf, inf};


    //declare vectors for the state and control at each node
    // std::vector<casadi::MX> X, U;
    for(int k=0; k<ns; ++k){
        // Local state
        // X.push_back( V.nz(casadi::Slice(offset,offset+nx)));
        if(k==0){
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        } else {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_init.begin(), x_init.end());

        // Local Control
        // U.push_back(V.nz(casadi::Slice(offset,offset+nu)));
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());
        v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    }

    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());

    // Bounds and initial guess
    casadi::Dict opts;

    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 200;
    opts["ipopt.linear_solver"] = "mumps";
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = 0;
    opts["ipopt.sb"] = "yes";

    Function solver = nlpsol("nlpsol", "ipopt", "double_pendulum.so", opts);

    std::map<std::string, casadi::DM> arg, res;
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;
    // arg["p"]  = traj;
    
    casadi::SX x_next = x + x_dot*time_step.as_seconds();
    casadi::Function F = casadi::Function("F",{x,u},{x_next},{"x","u"},{"x_next"});
    
    auto sim_time = mahi::util::seconds(5);
    double curr_sim_time = 0;
    
    std::vector<double> qA_opt,qB_opt,qA_dot_opt,qB_dot_opt;
    std::vector<double> TA_opt, TB_opt;
    
    mahi::util::Clock sim_clock;
    while (curr_sim_time < sim_time.as_seconds()){
        // std::cout << "running iteration " << iter++ << std::endl;
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
        arg["p"]  = traj;
        
        mahi::util::Clock my_clock;
        
        res = solver(arg);
        std::cout << my_clock.get_elapsed_time() << std::endl;
        
        std::vector<double> V_opt(res.at("x"));
        arg["x0"] = V_opt;
        
        std::vector<double> x_curr(V_opt.begin(),V_opt.begin()+nx);
        std::vector<double> u_curr(V_opt.begin()+(nx),V_opt.begin()+ (nx+nu));
        
        auto x_next = F(casadi::SXDict({{"x",x_curr},{"u",u_curr}}));

        std::vector<double> x_next_vec(x_next.at("x_next"));

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
    std::cout << "sim time: " << sim_clock.get_elapsed_time() << std::endl;
    // vector<double> V_opt(res.at("x"));
    
    // arg["x0"] = V_opt;


    // Extract the optimal state trajectory


    std::ofstream file;
    std::string filename = "new_multishoot_results.m";
    file.open(filename.c_str());
    file << "% Results from " __FILE__ << std::endl;
    file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
    file << std::endl;

    // Save results
    file << "t = linspace(0," << sim_time.as_seconds() << "," << (sim_time/time_step + 1) << ");" << std::endl;
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
