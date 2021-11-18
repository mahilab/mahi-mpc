#include <Mahi/Casadi/ModelGenerator.hpp>
#include <Mahi/Casadi/M.hpp>
#include <Mahi/Casadi/G.hpp>
#include <Mahi/Casadi/V.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using mahi::util::PI;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        // ("d,dll", "Generate new dll.")
        // ("s,s2", "use second timestep.")
        ("p,double_pendulum", "generate double pendulum instead of MOE")
        // ("m,ma27", "Runs using MA27 solver. Mumps otherwise")
        ("l,linear", "Generates linearized model.");
    
    auto result = options.parse(argc, argv);

    bool linear = result.count("linear") > 0;

    SX x, x_dot, u;
    std::string model_name;

    if (result.count("double_pendulum")){
        model_name = "double_pendulum";

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

        x = SX::vertcat({qA,qB,qA_dot,qB_dot});
        x_dot = SX::vertcat({qA_dot,qB_dot,qA_ddot,qB_ddot});
        
        // control vector
        u = SX::vertcat({TA,TB});
    }
    else{
        model_name = "moe";
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
        x = SX::vertcat({q0,q1,q2,q3,q0_dot,q1_dot,q2_dot,q3_dot});
        u = SX::vertcat({T0,T1,T2,T3});

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
        x_dot = vertcat(q_dot,q_d_dot);
    }
    if (linear) model_name = "linear_" + model_name;
    
    // Bounds on state
    std::vector<double> x_min(x.size1(),-inf);
    std::vector<double> x_max(x.size1(), inf);

    // Bounds for control
    std::vector<double> u_min(u.size1(),-inf);
    std::vector<double> u_max(u.size1(), inf);

    // settings for multiple shooting constructions
    mahi::util::Time time_step  = mahi::util::milliseconds(2);
    int num_shooting_nodes = 25;

    ModelParameters model_parameters(model_name, // name
                                     x.size1(),                // num_x
                                     u.size1(),                // num_u
                                     time_step,                // step_size
                                     num_shooting_nodes,       // num_shooting_nodes
                                     linear);                  // is_linear;                  

    // 
    ModelGenerator my_generator(model_parameters, x, x_dot, u);

    my_generator.create_model();
    my_generator.generate_c_code();
    my_generator.compile_model();

    return 0;
}
