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

    // set amplitude based on input option of double_pendulum
    double sin_amp;
    double sin_freq;
    if(result.count("double_pendulum")){
        sin_amp = 1.0;
        sin_freq = 1.0;
    } 
    else{
        sin_amp = 0.5;
        sin_freq = 0.25;
    }


    const int nx = model_control.model_parameters.num_x;

    mahi::util::Time sim_time = mahi::util::seconds(1.0);

    std::vector<double> state(model_control.model_parameters.num_x,0);
    std::vector<double> control(model_control.model_parameters.num_u,0);

    auto ext_x_dot_init = external(model_name + "_get_x_dot_init",model_name + "_linear_functions.so");

    std::vector<double> time_result;
    std::vector<std::vector<double>> x_result;
    std::vector<std::vector<double>> u_result;
    double avg_solve_time = 0;
    int cycle_counter = 0;
    
    // warm start
    double curr_traj_time = 0;
    for (auto i = 0; i < model_control.model_parameters.num_shooting_nodes; i++){
        for (auto j = 0; j < model_control.model_parameters.num_x; j++){
            // if it is one of the position states
            if (j < nx/2) traj.push_back( ((j % 2 == 0) ? 1.0 : -1.0) * sin_amp * sin(2*PI*sin_freq*curr_traj_time));
            // if it is one of the velocity states
            else          traj.push_back( (((j-nx/2) % 2 == 0) ? 1.0 : -1.0) * sin_amp * 2*PI*sin_freq*cos(2*PI*sin_freq*curr_traj_time));
        }
        curr_traj_time += model_control.model_parameters.step_size.as_seconds();
    }

    model_control.set_state(mahi::util::seconds(0), state, control, traj);
    model_control.start_calc();
    mahi::util::sleep(mahi::util::milliseconds(100));
    // end warm start
    mahi::util::Time sim_rate = mahi::util::microseconds(1000);
    mahi::util::Timer sim_clock(sim_rate);
    while (curr_sim_time < sim_time.as_seconds()){
        std::cout <<cycle_counter << "\n";
        // update trajectory with new desired
        mahi::util::Clock my_clock;
        curr_traj_time = curr_sim_time;
        traj.clear();
        for (auto i = 0; i < model_control.model_parameters.num_shooting_nodes; i++){
            for (auto j = 0; j < model_control.model_parameters.num_x; j++){
                // if it is one of the position states
                if (j < nx/2) traj.push_back( ((j % 2 == 0) ? 1.0 : -1.0) * sin_amp * sin(2*PI*sin_freq*curr_traj_time));
                // if it is one of the velocity states
                else          traj.push_back( (((j-nx/2) % 2 == 0) ? 1.0 : -1.0) * sin_amp * 2*PI*sin_freq*cos(2*PI*sin_freq*curr_traj_time));
            }
            curr_traj_time += model_control.model_parameters.step_size.as_seconds();
        }
    
        time_result.push_back(curr_sim_time);
        x_result.push_back(state);
        u_result.push_back(control);

        model_control.set_state(mahi::util::seconds(curr_sim_time), state, control, traj);
        auto control_result = model_control.control_at_time(mahi::util::seconds(curr_sim_time));
        control = control_result.u;
        // mahi::util::print("time: {}, state: {}, control: {}",control_result.time,state,control);

        std::vector<casadi::DM> lin_args = {state, control};
        std::vector<double> x_dot(ext_x_dot_init(lin_args)[0]);

        for (auto i = 0; i < state.size(); i++){
            state[i] += x_dot[i]*sim_rate.as_seconds();
        }

        auto solve_time = my_clock.get_elapsed_time().as_milliseconds();
        // std::cout << solve_time << ", ";
        avg_solve_time = (avg_solve_time*cycle_counter+(double)solve_time)/(cycle_counter+1);
        cycle_counter++;
        curr_sim_time = sim_clock.wait().as_seconds();
    }

    mahi::util::print("sim time: {:.2f} s, avg solve time: {:.2f} ms",sim_clock.get_elapsed_time().as_seconds(),avg_solve_time);

    model_control.stop_calc();

    // Extract the optimal state trajectory
    std::ofstream file;
    std::string filename = model_name + "_results_threaded.m";
    file.open(filename.c_str());

    // Save results
    mahi::util::Timestamp ts;
    file << "% " << ts.hh_mm_ss_mmm() << std::endl;
    // file << "t = linspace(0," << sim_time.as_seconds() << "," << x_result.size() << ");" << std::endl;
    file << "t = " << time_result << ";" << std::endl;
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
