#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <corecrt_math_defines.h>
#include <cmath>

//#include <boost/safe_numerics/safe_integer.hpp>

#include "file_h/complex.h"
#include "file_h/fraction.h"
#include "file_h/polynomial.h"
#include "file_h/matrix.h"

#include "file_h/counting_methods_1.h"
#include "file_h/counting_methods_2.h"
#include "file_h/Array_xy_To_.h"
#include "file_h/calc_computing_f.h"
#include "file_h/Spline.h"

template<typename T, typename Func>
std::vector<std::pair<T, T>> generatePointsLambda(int k, T x0, T step, Func F) {
    std::vector<std::pair<T, T>> points;
    for (int i = 0; i < k; ++i) {
        T x = x0 + i * step;
        points.emplace_back(x, F(x));
    }
    return points;
}

#define SHOW_INTERPOL_STAT(func,n,m, ...) \
    counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic(func,n+20,m, #func, __VA_ARGS__)

int main() {
    using namespace counting_methods_2::Polynomial_interpolation::nuton2;
#if 0
    //calc_computing_f::matrix_calc<double>();
//C:\Users\User\source\repos\cpp_rewriting_SPBGU\input_matrix.txt
    /*counting_methods::polinomial::polynomial_test();*/
    //counting_methods::nonlinear_system_with_simple_iterations::run_nnssi_with_setted_nonlinear_function();

    //counting_methods::nonlinear_system_with_the_tangent_method::nonlinsystem_tangent_method();
    //counting_methods::executeWithFileInput((counting_methods::gaus_method::gaus_solver_linear_sistem), "input_matrix.txt");
    //counting_methods::gaus_method::gaus_solver_linear_sistem();

    //counting_methods::nonlinear_system_with_simple_iterations::run_nnssi_with_setted_linear_function();

    counting_methods::holechi::example();
#endif
#define AU_LOG 7

#if AU_LOG==1
    
    std::vector<std::pair<int, int>> Array_xy = { {1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30} };
    Nuton<int> a(Array_xy);
    a.nuton_interpolation(6);
#elif AU_LOG==2
   
    auto TheTypeOfFunctionThatBeingInterpolated = [](double x) { return  1 + 10 * x + 10 * x * x + 10 * x * x * x + 10 * x * x * x * x;};

    auto Array_xy = generatePointsLambda(7, -4, 1, TheTypeOfFunctionThatBeingInterpolated);
    
    std::cout << nuton_interpolation(Array_xy)<<"\ndivided_difference0_to_k(Array_xy):"<< nuton_interpolation(Array_xy).get_deg() << "\n\n\n\n\n\n\n";

    /*auto Func2 = [](double x) { return -2479 - 1050 * x; };
    auto Array_xy_2 = generatePointsLambda(7, -4, 1, Func);

    std::cout << divided_difference0_to_k(Array_xy_2);*/


    
#elif AU_LOG==3

    int n = 10;
    double a = -2 * M_PI, b = 2 * M_PI;

    auto poly = N_optn(n, a, b, [](double x) { return std::cos(x) / std::sin(x) + x * x; });

    std::cout << "POLINOM:" << poly;
    std::cout << "\nMAX:" << RN_n(n, a, b, [](double x) { return std::cos(x) / std::sin(x) + x * x; });

#elif AU_LOG == 4
    //show_nuton();
    // 
    //int T0 = 1;
    //std::vector<double> array{1,2,3,4,5,6,7,8,9,10};

    //
   
    ////b.output_mode_set(output_mode::FULL);
    //polynomial<int> a(1);
    //for (uint64_t i = 0; i < array.size(); i++)
    //{
    //if (i!=T0){
    //    polynomial<int> b(1);
    //    b = 1; b = b >> 1; std::cout << b << '\n'; b = b + array[i]; std::cout << b << '\n';
    //    
    //    a = a * b;
    //    std::cout << a<<'\n';
    //    }
    //}
    //
    //
    //std::cout << a;

    std::vector<double> array={ 1,2,3,4,5,6 };
    std::cout<< "wk0" << w_k_T0(array, 1).output_mode_set(output_mode::ABBREVIATED);
    
    //static_cast<double(*)(double)>
    using T = double;
    //(w_k(array, 1))
    //draw_function<double(*)(double)>(std::cos);

    polynomial<double> pol = (w_k_T0(array, 1));
    std::vector<std::function<double(double)>> functions = {
        // Лямбда с захватом
        [k = 2.0](double x) { return k * x; },
        static_cast<double(*)(double)>(std::sin),
        pol
    };
    //draw_function([](float x) {return x * x + x * 2; },-25.0f,25.0f);
    draw_functions(functions);
    
    std::cout << pol(2);
#elif AU_LOG == 5

    std::vector<double> array = { 1,2,3,4,5,6 };
    polynomial<double> polynom = w_k_T0(array, -1);

    std::cout << "wk" << w_k_T0(array, -6).output_mode_set(output_mode::ABBREVIATED);
    std::cout << "\nwx0"<<(w_k_T0(array, 1))(1);

    std::vector<std::pair<double, double>> array_xy = {
        {0.5, polynom(0.5) },{2.5, polynom(2.5)},
        {4.5, polynom(4.5) },{6.5, polynom(6.5)},
        {8.5, polynom(8.5) },{10.5, polynom(10.5)}
    };


    auto Func = [](double x) { return  1 + 10 * x + 10 * x * x + 10 * x * x * x + 10 * x * x * x * x; };
    auto Array_xy = generatePoints_equally_sufficient_with_step_size(7, -4.0, 1.0, Func);

    std::cout <<"\nans:" << Lagrang_interpolation(Array_xy);
#elif AU_LOG == 6
    {
        polynomial<double>a; a.outm_E = output_mode::FULL;
        polynomial<double>b;
        std::vector <double> roots{ 1,3,0.5,0.6,12 };
        a.the_root_constructor(roots);
        auto V = (polynomialfunctions::plnm_roots(a.get_first_derrivate(), DBL_EPSILON));
        std::cout << a << "max|F(x)|=";
        std::cout << '\n';
        std::cout << (a) << "max|F(x)|=" << a.maximum_abs(-1111., +1111.) << '\n';
        b = nuton_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << '\n' << b << "max|F(x)|=" << b.maximum_abs(-1111., +1111.);
        b = Lagrang_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << '\n' << b << "max|F(x)|=" << b.maximum_abs(-1111., +1111.);
        b = Alternativ_nuton_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << '\n' << b << "max|F(x)|=" << b.maximum_abs(-1111., +1111.);
        b = Alternativ_Lagrang_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << '\n' << b << "max|F(x)|=" << b.maximum_abs(-1111., +1111.);

#define identifier (x-std::sin(x) - 0.25)
        //x-std::sin(x) - 0.25
        SHOW_INTERPOL_STAT(
            nuton_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            -M_PI / 4, M_PI / 4
        );
        SHOW_INTERPOL_STAT(
            Lagrang_interpolation<double>,
            10, 10,
            [](double x) { return identifier; },
            -M_PI / 4, M_PI / 4
        );
        SHOW_INTERPOL_STAT(
            Alternativ_nuton_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            -M_PI / 4, M_PI / 4
        );
        SHOW_INTERPOL_STAT(
            Alternativ_Lagrang_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            -M_PI / 4, M_PI / 4
        );



    }
    

#elif AU_LOG == 7

    std::vector<double> array = { 1,2,3,4,5,6 };
    polynomial<double> polynom = w_k_T0(array, -1);

    std::cout << "wk" << w_k_T0(array, -6).output_mode_set(output_mode::ABBREVIATED);
    std::cout << "\nwx0" << (w_k_T0(array, 1))(1);

    std::vector<std::pair<double, double>> array_xy = {
        {0.5, polynom(0.5) },{2.5, polynom(2.5)},
        {4.5, polynom(4.5) },{6.5, polynom(6.5)},
        {8.5, polynom(8.5) },{10.5, polynom(10.5)}
    };
    
    Spline<double> a(-5,5,10);
    std::cout << a;
    Spline_build<2,1>(array_xy);


    
#endif
    
    system("pause");
    return 0;
}
