#include <Mahi/Casadi/ModelParameters.hpp>
#include <casadi/casadi.hpp>

ModelParameters::ModelParameters(std::string name_, int num_x_, int num_u_, mahi::util::Time step_size_, size_t num_shooting_nodes_, bool is_linear_, std::vector<double> u_min_, std::vector<double> u_max_, std::vector<double> x_min_,  std::vector<double> x_max_):
    name(name_),
    num_x(num_x_),
    num_u(num_u_),
    step_size(step_size_),
    num_shooting_nodes(num_shooting_nodes_),
    is_linear(is_linear_),
    x_min(x_min_),
    u_min(u_min_),
    x_max(x_max_),
    u_max(u_max_)
    {
        timespan = mahi::util::microseconds(step_size.as_microseconds()*num_shooting_nodes);
        
        if (x_min.empty()) x_min = std::vector<double>(num_x,-10e30);
        if (x_max.empty()) x_max = std::vector<double>(num_x,10e30);
        if (u_min.empty()) u_min = std::vector<double>(num_u,-10e30);
        if (u_max.empty()) u_max = std::vector<double>(num_u,10e30);
    }

void to_json(mahi::util::json& j, const ModelParameters& p) {
    j = mahi::util::json{{"name", p.name}, 
                         {"timespan", p.timespan.as_microseconds()}, 
                         {"step_size", p.step_size.as_microseconds()}, 
                         {"num_x", p.num_x},
                         {"num_u", p.num_u},
                         {"num_shooting_nodes", p.num_shooting_nodes},
                         {"x_min", p.x_min},
                         {"u_min", p.u_min},
                         {"x_max", p.x_max},
                         {"u_max", p.u_max},
                         {"dll_filepath", p.dll_filepath},
                         {"is_linear", p.is_linear}};
}

void from_json(const mahi::util::json& j, ModelParameters& p) {
    j.at("name").get_to(p.name);
    p.timespan = mahi::util::microseconds(j.at("timespan"));
    p.step_size = mahi::util::microseconds(j.at("step_size"));
    j.at("num_x").get_to(p.num_x);
    j.at("num_u").get_to(p.num_u);
    j.at("num_shooting_nodes").get_to(p.num_shooting_nodes);
    j.at("x_min").get_to(p.x_min);
    j.at("u_min").get_to(p.u_min);
    j.at("x_max").get_to(p.x_max);
    j.at("u_max").get_to(p.u_max);
    j.at("dll_filepath").get_to(p.dll_filepath);
    j.at("is_linear").get_to(p.is_linear);

    for (size_t i = 0; i < p.x_min.size(); i++){
        if (p.x_min[i] == -10e30) p.x_min[i] = -casadi::inf;
        if (p.x_max[i] ==  10e30) p.x_max[i] =  casadi::inf;
    }
    
    // std::cout << p.name;
}