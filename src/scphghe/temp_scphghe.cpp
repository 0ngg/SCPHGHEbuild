#include"temp_struct_cfd.h"
#include<conio.h>

int main(int argc, char** argv)
{
    // args
    // 0 mesh name, 1 output name, 2 source file, 3 prop_solid file
    // 4 double P_init, 5 double T_init, 6 double W_init
    // 7 double step_length, 8 double under_relax, 9 double min_residual, 10 int max_iter;
    std::cout << "Initializing..." << std::endl;
    std::string mesh_arg(argv[0]); std::string output_arg(argv[1]); std::string source_arg(argv[2]); std::string prop_solid_arg(argv[3]);
    double p_arg = std::stod(argv[4]); double T_arg = std::stod(argv[5]); double W_arg = std::stod(argv[6]);
    double step_length_arg = std::stod(argv[7]); double under_relax_arg = std::stod(argv[8]); double min_res_arg = std::stod(argv[9]);
    int max_iter_arg = std::stoi(argv[10]);
    cfd::user* main_user(new cfd::user(p_arg, T_arg, W_arg, source_arg, prop_solid_arg));
    std::cout << "User obj. successfully built." << std::endl;
    cfd::scheme* main_scheme(new cfd::scheme(mesh_arg, main_user));
    std::cout << "Scheme obj. successfully built." << std::endl;
    _getch();
    return 0;
};