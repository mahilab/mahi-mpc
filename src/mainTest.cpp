// #include <casadi/casadi.hpp>
#include <casadi/casadi.hpp>
#include <iostream>
#include <fstream>
#include <ctime>
// #include <Mahi/Util/Print.hpp>
// #include <Mahi/Util/Math/Constants.hpp>
#include <Mahi/Util.hpp>

int main() {
    casadi::Slice all;
    int T = 10;
    int N = 200;
    double L = 1.0;
    double m = 1.0;
    double g = 9.81;
    double dt = 10.0/200.0;

    // mahi::util::print("Test");
    //casadi::MX x =  casadi::MX::sym("x",2,2);
    casadi::MX qA = casadi::MX::sym("qA",1,1);
    casadi::MX qB = casadi::MX::sym("qB",1,1);
    casadi::MX qA_dot = casadi::MX::sym("qA_dot",1,1);
    casadi::MX qB_dot = casadi::MX::sym("qB_dot",1,1);

    casadi::MX x = casadi::MX::vertcat({qA,qB,qA_dot,qB_dot});

    casadi::MX TA = casadi::MX::sym("TA",1,1);
    casadi::MX TB = casadi::MX::sym("TB",1,1);
    casadi::MX u = casadi::MX::vertcat({TA,TB});
    //double test = pow(2,2);
    casadi::MX qA_ddot =  -(TA - TB - TB*cos(qB) + pow(L,2)*m*pow(qA_dot,2)*sin(qB) + pow(L,2)*m*pow(qB_dot,2)*sin(qB) - 2*L*g*m*cos(qA) + pow(L,2)*m*pow(qA_dot,2)*cos(qB)*sin(qB) + 2*pow(L,2)*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(pow(L,2)*m*(pow(cos(qB),2) - 2));
    casadi::MX qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*pow(L,2)*m*pow(qA_dot,2)*sin(qB) + pow(L,2)*m*pow(qB_dot,2)*sin(qB) - 2*L*g*m*cos(qA) + 2*pow(L,2)*m*pow(qA_dot,2)*cos(qB)*sin(qB) + pow(L,2)*m*pow(qB_dot,2)*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*pow(L,2)*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*pow(L,2)*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(pow(L,2)*m*(pow(cos(qB),2) - 2));
    casadi::MX ode = casadi::MX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});
    
    casadi::Function f = casadi::Function("f",{x,u},{ode},{"x","u"},{"ode"});
    

    casadi::Dict intg_options;

    intg_options["tf"] = 10.0/200.0;
    intg_options["simplify"] = true;
    intg_options["number_of_finite_elements"] = 4;

    

    casadi::MXDict dae;
    dae["x"] = x;
    dae["p"] = u;
    auto fRes =\
    f({{"x",x},{"u",u}});
    auto fTest = fRes["ode"];
    dae["ode"] = fTest;//(x,u);

    // casadi::Function intg = casadi::integrator("intg","rk",f,intg_options);
    casadi::Function intg = casadi::integrator("intg","rk",dae,intg_options);
    
    auto res =\
    intg({{"x0",x},{"p",u}});

    auto x_next = res["xf"];


    casadi::Function F = casadi::Function("F",{x,u},{x_next},{"x","u"},{"x_next"});
    auto rkRes =\
    F({{"x",x},{"u",u}});
    auto myNext = rkRes["x_next"];
    //std::cout << myNext.rows() << std::endl;
    


    casadi::Opti opti;
    auto x_opti = opti.variable(4,N+1);
    auto u_opti = opti.variable(2,N);
    auto Q = opti.parameter(4,4);
    opti.set_value(Q,casadi::DM::eye(4));
    
    casadi::MX R = opti.parameter(2,2);
    opti.set_value(R,casadi::DM::eye(2));

    casadi::MX x0param = opti.parameter(4,1);
    opti.set_value(x0param,casadi::DM::zeros(4,1));
    // casadi::MX t0 = opti.parameter(1,1);

    casadi::MX cost = 0;
    double desired_pos = 0.0;
    double desired_vel = 0.0;
    double current_t = 0.0;
    std::vector<double> desired_state(4,0);
    casadi::MX error = casadi::MX::sym("error",(4,1));
    for (size_t i = 0; i < N; i++)
    {
        current_t = i*dt;
        desired_pos = mahi::util::PI/2*sin(current_t);
        desired_vel = mahi::util::PI/2*cos(current_t);
        desired_state = {desired_pos,desired_pos,desired_vel,desired_vel};
        error = x_opti(all,i)-desired_state;
        cost = cost + mtimes(error.T(),mtimes(Q,error));// + mtimes(u_opti(all,i).T(),mtimes(R,u_opti(all,i)));
    }
    desired_pos = mahi::util::PI/2*sin(N*dt);
    desired_vel = mahi::util::PI/2*cos(N*dt);
    desired_state = {desired_pos,desired_pos,desired_vel,desired_vel};
    error = x_opti(all,N)-desired_state;
    cost = cost + mtimes(error.T(),mtimes(Q,error));

    casadi::Function J = casadi::Function("J",{x_opti,u_opti},{cost},{"x_opti","u_opti"},{"cost"});

    auto resCost =\
    J({{"x_opti",x_opti},{"u_opti",u_opti}});
    auto myCost = resCost["cost"];

    opti.minimize(myCost);

    casadi::MXDict tempVar;
        
    for (size_t i = 0; i < N; i++)
    {
        tempVar =\
         F({{"x",x_opti(all,i)},{"u",u_opti(all,i)}});
        opti.subject_to(x_opti(all,i+1)==tempVar["x_next"]);
    }

    opti.subject_to(x_opti(all,0) == x0param);

    casadi::Dict mySolverOpts;
    // mySolverOpts["verbose"] = false;
    // mySolverOpts["expand"] = true;
    // mySolverOpts["ipopt.print_level"] = 5;
    // mySolverOpts["ipopt.max_iter"] = 200;
    // mySolverOpts["error_on_fail"] = false;
    // //mySolverOpts["print_level"] = 1;
    // opti.solver("ipopt",mySolverOpts,casadi::Dict());//mySolverOpts)

    mySolverOpts["expand"] = true;
    mySolverOpts["qpsol"] = "nlpsol";
    mySolverOpts["qpsol_options.nlpsol"] = "ipopt";
    mySolverOpts["qpsol_options.error_on_fail"] = false;
    mySolverOpts["qpsol_options.nlpsol_options.ipopt.max_iter"] = 2;
    mySolverOpts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0;
    mySolverOpts["qpsol_options.nlpsol_options.print_time"] = 0;
    opti.solver("sqpmethod",mySolverOpts);
    // casadi::Function M = opti.to_function("M",{x0param},{u_opti},{"x0"},{"u"});
    // M.generate_dependencies("nlp.c");

    // std::cout << "Display please" << std::endl;   

    mahi::util::Clock clock;
    
    auto sol = opti.solve();
    mahi::util::Time solveTime = clock.get_elapsed_time();
    mahi::util::print("{}",solveTime);

    std::ofstream file;
    std::string filename = "casadi_test_results.m";
    file.open(filename.c_str());
    file << "% Results file from " __FILE__ << std::endl;
    file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
    file << std::endl;

    // Save results to file
    file << "t = linspace(0," << T << "," << N << "+1);" << std::endl;
    file << "qA = " << std::vector<double>(sol.value(x_opti(0,all))) << ";" << std::endl;
    file << "qB = " << std::vector<double>(sol.value(x_opti(1,all))) << ";" << std::endl; 
    file << "qA_dot = " << std::vector<double>(sol.value(x_opti(2,all))) << ";" << std::endl;
    file << "qB_dot = " << std::vector<double>(sol.value(x_opti(3,all))) << ";" << std::endl;
    file << "T_A = " << std::vector<double>(sol.value(u_opti(0,all))) << ";" << std::endl;
    file << "T_B = " << std::vector<double>(sol.value(u_opti(1,all))) << ";" << std::endl;

    file << "figure;" << std::endl;
    file << "hold on;" << std::endl;
    file << "plot(t,qA);" << std::endl;
    file << "plot(t,qB);" << std::endl;
    file << "xlabel('Time (s)');" <<std::endl;
    file << "ylabel('Position (rad)');" << std::endl;
    file << "legend('qA','qB');" << std::endl;

    
    return 0;
}