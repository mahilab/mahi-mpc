#include <Mahi/Mpc/ModelGenerator.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using mahi::util::PI;
using namespace mahi::mpc;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("l,linear", "Generates linearized model.");
    
    auto result = options.parse(argc, argv);

    bool linear = result.count("linear") > 0;

    SX x, x_dot, u;
    std::string model_name;

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
    if (linear) model_name = "linear_" + model_name;
    else model_name = "nonlinear_" + model_name;
    
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
