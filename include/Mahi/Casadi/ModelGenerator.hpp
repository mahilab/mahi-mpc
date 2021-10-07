#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>

class ModelGenerator
{
private:
    std::string m_name; // name of the model for useful outputs
    mahi::util::Time m_timespan;
    mahi::util::Time m_step_size;
    casadi::SX m_x; // number of control inputs
    casadi::SX m_x_dot; // number of control inputs
    casadi::SX m_u; // number of control inputs
    casadi::SX m_quad; // number of control inputs
    size_t m_num_x; // number of states
    size_t m_num_u; // number of control inputs
    size_t m_num_shooting_nodes; // number of shooting_ndes
    std::vector<double> m_x_min; // minimum of states
    std::vector<double> m_u_min; // minimum of control inputs
    std::vector<double> m_x_max; // maximum of states
    std::vector<double> m_u_max; // maximum of control inputs
    std::string m_dll_filepath;
    std::string m_c_file_filepath;
    bool m_linear;

    casadi::Function m_solver;
    
public:
    
    ModelGenerator(std::string name, casadi::SX x, casadi::SX x_dot, casadi::SX u, mahi::util::Time final_time, mahi::util::Time step_size, std::vector<double> u_min, std::vector<double> u_max, std::vector<double> x_min, std::vector<double> x_max, bool linear = false);
    ~ModelGenerator();
    void create_model();
    void generate_c_code(std::string filepath);
    void compile_model(std::string filepath);

};