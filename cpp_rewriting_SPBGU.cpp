#include "file_h\project_libraries.h"



int main() {
    matrix< double> T{ {1,5,9}};
    std::cout << "T :\n" << T<<"\n";
    //T=T.transpose();
    std::cout << "coeff :\n" << counting_methods_3::СoefficNewtonCotes(T,static_cast<double>(1), static_cast<double>(9)).set_output_mode(output_mode::FULL) << "\n";
    //std::cout << "T :\n" << T << "\n";
    auto F = [](double x) {return x * x-1; };
    std::cout << "F(3)=" << F(3) << '\n';
    GetaPartialIntegralSum_Singlethreaded([](double x) {return x*x*4+2*x*x; }, -1, 1, 1, { {0,0.25,0.5,0.75,1} });
    system("pause");
    return 0;
}
