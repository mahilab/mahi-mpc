#include <casadi/casadi.hpp>
#include <Mahi/Util.hpp>

using namespace std;
using namespace casadi;

int main(int argc, char const *argv[])
{
    map<string, DM> arg, res;

    // Bounds and initial guess
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;

    
    
    // Solve the problem
    res = solver(arg);
    // The optimal solution
    vector<double> V_opt(res.at("x"));
    
    // Extract the optimal state trajectory
    vector<double> qA_opt(ns+1),qB_opt(ns+1),qA_dot_opt(ns+1),qB_dot_opt(ns+1);

    for(int i=0; i<=ns; ++i){
        qA_opt[i] = V_opt.at(i*(nx+2));
        qB_opt[i] = V_opt.at(1+i*(nx+2));
        qA_dot_opt[i] = V_opt.at(2+i*(nx+2));
        qB_dot_opt[i] = V_opt.at(3+i*(nx+2));
    }

    // Get the optimal control
    vector<double> TA_opt(ns), TB_opt(ns);
    for(int i=0; i<ns; ++i){
        TA_opt[i] = V_opt.at(nx + i*(nx+2));
        TB_opt[i] = V_opt.at(nx + 1 + i*(nx+2));
    }

    ofstream file;
    string filename = "my_multishoot_results.m";
    file.open(filename.c_str());
    file << "% Results from " __FILE__ << endl;
    file << "% Generated " __DATE__ " at " __TIME__ << endl;
    file << endl;

    // Save results
    file << "t = linspace(0,10," << ns << "+1);" << endl;
    file << "qA = " << qA_opt << ";" << endl;
    file << "qB = " << qB_opt << ";" << endl;
    file << "qA_dot = " << qA_dot_opt << ";" << endl;
    file << "qB_dot = " << qB_dot_opt << ";" << endl;
    file << "TA = " << TA_opt << ";" << endl;
    file << "TB = " << TB_opt << ";" << endl;

    file << "figure;" << endl;
    file << "hold on;" << endl;
    file << "plot(t,qA);" << endl;
    file << "plot(t,qB);" << endl;
    file << "xlabel('Time (s)');" << endl;
    file << "ylabel('Position (rad)');" << endl;
    file << "legend('qA','qB');" << endl; 
    cout << "finished" << endl;
    return 0;
}
