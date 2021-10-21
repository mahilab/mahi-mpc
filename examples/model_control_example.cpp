#include <Mahi/Util.hpp>
#include <Mahi/Casadi/ModelControl.hpp>

using mahi::util::PI;
using namespace casadi;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("p,double_pendulum", "generate double pendulum instead of MOE")
        ("l,linear", "Generates linearized model.");

    auto result = options.parse(argc, argv);

    std::string model_name = result.count("double_pendulum") ? "double_pendulum" : "moe";
    if (result.count("linear")) model_name = "linear_" + model_name;
    std::cout << "Loading " << model_name << std::endl;
    casadi::Dict solver_opts;
    ModelControl model_control(model_name, solver_opts);

    std::vector<double> traj;

    double curr_sim_time = 0.0;

    const int nx = model_control.model_parameters.num_x;

    mahi::util::Time sim_time = mahi::util::seconds(1.0);

    std::vector<double> state(model_control.model_parameters.num_x,0);
    std::vector<double> control(model_control.model_parameters.num_u,0);

    auto ext_x_dot_init = external(model_name + "_get_x_dot_init",model_name + "_linear_functions.so");

    std::vector<std::vector<double>> x_result;
    std::vector<std::vector<double>> u_result;
    double avg_solve_time = 0;
    int cycle_counter = 0;
    mahi::util::Clock sim_clock;
    while (curr_sim_time < sim_time.as_seconds()){
        // update trajectory with new desired
        mahi::util::Clock my_clock;
        double curr_traj_time = curr_sim_time;
        traj.clear();
        for (auto i = 0; i < model_control.model_parameters.num_shooting_nodes; i++){
            for (auto j = 0; j < model_control.model_parameters.num_x; j++){
                // if it is one of the position states
                if (j < nx/2) traj.push_back( ((j % 2 == 0) ? 1.0 : -1.0) * sin(2*PI*curr_traj_time));
                // if it is one of the velocity states
                else          traj.push_back( (((j-nx/2) % 2 == 0) ? 1.0 : -1.0) * 2*PI*cos(2*PI*curr_traj_time));
            }
            curr_traj_time += model_control.model_parameters.step_size.as_seconds();
        }

        x_result.push_back(state);
        u_result.push_back(control);

        if (cycle_counter % 1 == 0){
            auto result = model_control.calc_u(mahi::util::seconds(curr_sim_time), state, control, traj);
        }
        auto control_result = model_control.control_at_time(mahi::util::seconds(curr_sim_time));

        control = control_result.u;

        std::vector<casadi::DM> lin_args = {state, control};
        std::vector<double> x_dot(ext_x_dot_init(lin_args)[0]);

        for (auto i = 0; i < state.size(); i++){
            state[i] += x_dot[i]*model_control.model_parameters.step_size.as_seconds();
        }

        curr_sim_time += model_control.model_parameters.step_size.as_seconds();
        auto solve_time = my_clock.get_elapsed_time().as_milliseconds();
        // std::cout << solve_time << ", ";
        avg_solve_time = (avg_solve_time*cycle_counter+(double)solve_time)/(cycle_counter+1);
        cycle_counter++;
    }

    mahi::util::print("sim time: {:.2f} s, avg solve time: {:.2f} ms",sim_clock.get_elapsed_time().as_seconds(),avg_solve_time);

    // Extract the optimal state trajectory
    std::ofstream file;
    std::string filename = model_name + "_results.m";
    file.open(filename.c_str());

    // Save results
    file << "t = linspace(0," << sim_time.as_seconds() << "," << x_result.size() << ");" << std::endl;
    for (size_t i = 0; i < x_result[0].size(); i++){
        file << "q" << i << " = [";
        for (size_t j = 0; j < x_result.size(); j++){
            file << x_result[j][i];
            if (j != x_result.size()-1) file << ", ";
        }
        file <<  "];" << std::endl;
    }

    for (size_t i = 0; i < u_result[0].size(); i++){
        file << "T" << i << " = [";
        for (size_t j = 0; j < u_result.size(); j++){
            file << u_result[j][i];
            if (j != u_result.size()-1) file << ", ";
        }
        file <<  "];" << std::endl;
    }

    file << "figure;" << std::endl;
    file << "hold on;" << std::endl;
    for (size_t i = 0; i < x_result[0].size(); i++){
        file << "plot(t,q" << i << ");" << std::endl;
    }
    file << "xlabel('Time (s)');" << std::endl;
    file << "ylabel('Position (rad)');" << std::endl;
    file << "legend(";
    for (size_t i = 0; i < x_result[0].size(); i++){
        file << "'q" << i;
        if (i != x_result[0].size()-1) file << "', ";
    }
    file << "');" << std::endl; 

    file << "figure;" << std::endl;
    file << "hold on;" << std::endl;
    for (size_t i = 0; i < u_result[0].size(); i++){
        file << "plot(t,T" << i << ");" << std::endl;
    }
    file << "xlabel('Time (s)');" << std::endl;
    file << "ylabel('Position (rad)');" << std::endl;
    file << "legend(";
    for (size_t i = 0; i < u_result[0].size(); i++){
        file << "'T" << i;
        if (i != u_result[0].size()-1) file << "', ";
    }
    file << "');" << std::endl; 
    
    return 0;
}
