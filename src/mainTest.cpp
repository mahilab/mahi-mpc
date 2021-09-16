// #include <casadi/casadi.hpp>
#include <C:/MatlabHelpers/Casadi/casadi-windows-matlabR2016a-v3.5.5/include/casadi/casadi.hpp>
#include <iostream>
#include <fstream>
#include <Mahi/Util/Print.hpp>

int main() {
    casadi::Slice all;
    int T = 10;
    int N = 200;
    double L = 1.0;
    double m = 1.0;
    double g = 9.81;

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
    //auto Q = opti.parameter(4,4);
    //opti.set_value(Q,casadi::DM::eye(4));

    

    //casadi::MX R = opti.parameter(2,2);
    //opti.set_value(R,casadi::DM::eye(2));

    casadi::MX x0param = opti.parameter(4,1);
    opti.set_value(x0param,casadi::DM::zeros(4,1));
    // casadi::MX t0 = opti.parameter(1,1);

    


    // double tCurrent;
    // for (size_t i = 0; i < N; i++)
    // {
    //     /* code */
    //     cost = cost + 1000*((x(i).T())*Q*x(i)) + u(i).T()*R*u(i);
    // }
    // cost = cost + x(N).T()*Q*x(N);

    auto cost = casadi::MX::sumsqr(x_opti)+ casadi::MX::sumsqr(u_opti);
    casadi::Function J = casadi::Function("J",{x_opti,u_opti},{cost},{"x_opti","u_opti"},{"cost"});




    auto resCost =\
    J({{"x_opti",x_opti},{"u_opti",u_opti}});
    auto myCost = resCost["cost"];


    opti.minimize(myCost);

    casadi::MXDict tempVar;
        
    for (size_t i = 0; i < N; i++)
    {
        /* code */
        // std::cout << "Count this" << std::endl;
        // casadi::MX ode = casadi::MX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});
        tempVar =\
         F({{"x",x_opti(all,i)},{"u",u_opti(all,i)}});
        // auto tempVar2 = tempVar["x_next"];
        opti.subject_to(x_opti(all,i+1)==tempVar["x_next"]);
    }



    opti.subject_to(x_opti(all,0) == x0param);

    //std::cout << "Display please" << std::endl;   


    opti.solver("ipopt");

    //std::cout << "Display please" << std::endl;   


    auto sol = opti.solve();
    
    std::cout << "Display please" << std::endl;   

    return 0;
}