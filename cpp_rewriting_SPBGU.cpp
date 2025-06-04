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
#include "file_h/eigenvalues.h"

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
    counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic(func,n,m+50, #func, __VA_ARGS__)


//

void test_solve_system() {
    // Создаем диагональную матрицу 3x3
    matrix<double> A = matrix<double>::zeros(3, 3);
    A[0][0] = 2.0;
    A[1][1] = 3.0;
    A[2][2] = 4.0;

    // Вектор правой части
    matrix<double> b = matrix<double>::ones(3, 1);
    b[0][0] = 4.0;
    b[1][0] = 6.0;
    b[2][0] = 8.0;

    // Решаем систему: A * x = b
    matrix<double> x = matrixfunction::solve_system(A, b);

    // Ожидаемое решение: [2.0, 2.0, 2.0]
    std::cout << "Solution:\n" << x << std::endl;
}

void test_solve_system_complex() {
    // Матрица A (недиагональная, требует перестановок)
    matrix<double> A = matrix<double>::zeros(3, 3);
    A[0][0] = 0.0;  A[0][1] = 2.0;  A[0][2] = 1.0;
    A[1][0] = 1.0;  A[1][1] = 1.0;  A[1][2] = 1.0;
    A[2][0] = 2.0;  A[2][1] = 0.0;  A[2][2] = 3.0;
    //std::cout << A.set_output_mode(output_mode::ABBREVIATED) << "\n";




    // Вектор правой части: b = [5, 6, 13]
    matrix<double> b = matrix<double>::ones(3, 1);
    b[0][0] = 5.0;
    b[1][0] = 6.0;
    b[2][0] = 13.0;

    // Решаем систему: A * x = b
    matrix<double> x = matrixfunction::solve_system(A, b);


    // Ожидаемое решение: x = [1, 2, 3]
    std::cout << "Computed solution:\n" << x << std::endl;

    // Проверка LU-разложения: L * U = P * A
    auto lup = A.LUP();
    b.set_output_mode(output_mode::ABBREVIATED);
    auto T = (lup.P * b);T.set_output_mode(output_mode::ABBREVIATED);
    std::cout << "b:\n" << T.set_output_mode(output_mode::ABBREVIATED) << "\n";
    matrix<double> LU = (lup.L * lup.U);LU.set_output_mode(output_mode::ABBREVIATED);
    /*std::cout << "l:\n" << ((LU.get_output_mode()) == output_mode::FULL) << "\n";*/
    matrix<double> PA = (lup.P * A);PA.set_output_mode(output_mode::ABBREVIATED);
    std::cout << "L * U:\n" << LU << "\nP * A:\n" << PA << std::endl;

}

void Z5_2(int Size=4) {
    std::cout << "Z5_2:\n\n";
    matrix<double> A;

    A = matrix<double>::randomDiagonal(Size, -100, 100);
    std::cout << "A:\n" << A << "\n";
    matrix<double> T = (matrix<int>::random(Size, Size, -100, 100));
    A = T * A * (T.inverse_M());
    std::cout << "T * A *(T^-1):\n" << A << "\n";

    std::vector<double> shifts = { 0.9, -1.5, 4.5 ,20 };
    std::vector<std::pair<double, matrix<double>>> result = matrixfunction::inverse_power_method_with_shifts(A, shifts);
    for (auto p : result)
    {
        std::cout << "Eigenvalue: " << p.first << "\nEigenvector:\n" << p.second << "\n";
    }
    matrix<double> vec0 = matrix<double>::ones(Size, 1);
    // std::cout << vec0 << "\n";
    auto result_ = matrixfunction::inverse_power_method_with_shift(A, 1.5, vec0);
    std::cout << "Eigenvalue: " << result_.first << "\nEigenvector:\n" << result_.second;

}
void Z5_1_2_const() {
    std::cout << "Z5_2_const:\n\n";
    matrix<double> A=matrix<double>::zeros(4,4);
    A[0][0] = 120; 
     A[1][1] = 12.0;  
     A[2][2] = 1.2; 
     A[3][3] = 0.012;
    std::cout << "A:\n" << A << "\n";
    
    matrix<double> T(4, 4);
    T[0][0] = 76;  T[0][1] = -56;  T[0][2] = -44;  T[0][3] = -80;
    T[1][0] = 49;  T[1][1] = 49;  T[1][2] = 100;  T[1][3] = 44;
    T[2][0] = 85;  T[2][1] = -37;  T[2][2] = 93;  T[2][3] = 20;
    T[3][0] = 39;  T[3][1] = 95;  T[3][2] = -2;  T[3][3] = -22;

    std::cout << "T :\n" << T << "\n";
    A = T * A * (T.inverse_M());
    std::cout << "T * A *(T^-1):\n" << A << "\n";

    
    std::vector<double> shifts = { 1000, 10, 1 ,0 }; 
    std::vector<std::pair<double, matrix<double>>> result = matrixfunction::inverse_power_method_with_shifts(A, shifts);
    std::cout << "inverse_with_shifts_Power_meth:\n\n";
    for (auto p : result)
    {
        std::cout << "Eigenvalue: " << p.first << "\nEigenvector:\n" << p.second << "\n";
    }
    
    std::cout << "matrixfunction::power_method_average:\n\n";
    auto Ans = matrixfunction::power_method_average(A, matrix<double>::ones(A.getcol(), 1), 1e-10, 10000);
    //matrixfunction::power_method_max(A, matrix<double>::ones(A.getcol(), 1),1e-10,10000);
    std::cout << "Eigenvalue: " << Ans.first << "\nEigenvector:\n" << matrixfunction::simplify_eigenvector(Ans.second) << "\n";
}

void Z5_plus(int Size = 4) {
    std::cout << "Z5_plus:\n\n";
    //fraction<polynomial<double>>
    matrix<polynomial<double>> A;

    A = matrix<double>::randomDiagonal(Size, -100, 100);
    std::cout << "A:\n" << A << "\n";

    matrix<polynomial<double>> T = (matrix<int>::random(Size, Size, -100, 100));
    std::cout << "T:\n" << T << "\n";
    A = T * A * (T.inverse_M());
    std::cout << "T * A *(T^-1):\n" << A << " \n";
    //std::cout<< std::fixed << A.determinant() << "\n";

    matrix<polynomial<double>> A_ = A - matrix<polynomial<double>>::eye(Size) * (polynomial<double>(1) * (polynomial<double>(1) >> 1));
    
    //std::cout << "A_:\n" << A_ << "\n";
    uint64_t j = 0;
    for (auto& i: A_.determinant().plnm_roots())
    {
        j++;
        std::cout << "root_["<<j<<"]= "<< i.first << "\n";
    }
    



}

void Z5_qr(int Size = 4) {

    std::cout << "\n\nZ5_qr:\n\n";
    //fraction<polynomial<double>>
    matrix<double> A;

    A = matrix<double>::randomDiagonal(Size, -100, 100);
    std::cout << "A:\n" << A << "\n";

    matrix<double> T = (matrix<int>::random(Size, Size, -100, 100));
    std::cout << "T:\n" << T << "\n";
    A = T * A * (T.inverse_M());
    std::cout << "T * A *(T^-1):\n" << A << " \n  ";
    //std::cout << "A_:\n" << A_ << "\n";
    uint64_t j = 0;
    for (auto& i : matrixfunction::compute_eigenvalues_3_qr(A))
    {
        j++;
        std::cout << "eig[" << j << "]= "  << i<< "\n";
    }


}
//
// 
// 
template<typename T>
std::vector<std::pair<T, T>> addNoiseToPoints(
    const std::vector<std::pair<T, T>>& points,
    int measurements_per_point = 3,
    T noise_level = T(0.1)
) {


    std::vector<std::pair<T, T>> noisy_points;
    noisy_points.reserve(points.size() * measurements_per_point);

    std::random_device rd;
    std::mt19937 gen(rd());

    for (const auto& point : points) {
        const T x = point.first;
        const T y_clean = point.second;

        std::uniform_real_distribution<T> dist(-noise_level, noise_level);

        for (int i = 0; i < measurements_per_point; ++i) {
            T noise = dist(gen);
            T y_noisy = y_clean + noise;

            noisy_points.emplace_back(x, y_noisy);
        }
    }

    return noisy_points;
}
void Z6_1(int size,int degree_of_the_polynomial) {
#define loc_identifier1 x*x + 20*x*x*x+ 444*x*x*x*x*x 
#define loc_identifier2 (x-std::sin(x) - 0.25)
#define identifier loc_identifier1
#define the_left_border -M_PI / 4
#define the_right_border M_PI / 4
//#define SIZE_ 20
//#define degree_of_the_polynomial 7


    std::vector<std::pair<double, double>> clean_points = counting_methods_2::Polynomial_interpolation::nuton2::generatePoints_equally_sufficient_(size, the_left_border, the_right_border, [](double x) { return identifier; });
    auto noisy_points = addNoiseToPoints(clean_points, 3, 0.2);
    noisy_points.insert(noisy_points.end(), clean_points.begin(), clean_points.end());
    std::sort(noisy_points.begin(), noisy_points.end());
    for (auto& I : noisy_points) {
        std::cout << "<" << I.first << ";" << I.second << ">\n";
    }
    std::cout << "\npolynomial function :\n"<< "x*x + 20*x*x*x+ 444*x*x*x*x*x\n";

    std::cout << "\nNormal equations polynomial:\n" << counting_methods_2::aproximate::least_squares_normal(noisy_points, degree_of_the_polynomial);
    std::cout << "\nOrthogonal polynomial:\n" << counting_methods_2::aproximate::least_squares_orthogonal(noisy_points, degree_of_the_polynomial);
}
//
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



#define AU_LOG 0
    Z5_qr();

    // // // //
#if AU_LOG==0
    Z6_1(20, 7);

#define identifier (x-std::sin(x) - 0.25)
#define the_left_border -M_PI / 4
#define the_right_border M_PI / 4
    using SplineInterpolatorFunc = Spline<double>(*)(std::vector<std::pair<double, double>>);

    //x-std::sin(x) - 0.25
    counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic<double, SplineInterpolatorFunc>(
        Spline_interpolator<3, 2, double>, // interpolator
        50, 50,                      // n, m_
        "Spline_interpolator<3, 2, double>",
        [](double x) { return  identifier; },
        the_left_border, the_right_border,
        5
    );


    counting_methods_2::aproximate::show_aproximate_statistic([](double x) { return (x - std::sin(x) - 0.25); });

// // // // // // // // //


#elif AU_LOG==1

#if 0
    Z5_2(4);
    Z5_1_2_const();
    Z5_plus();
    Z5_qr();


    //part 3




#else
    Z6_1(20,7);


    counting_methods_2::aproximate::show_aproximate_statistic([](double x) { return (x - std::sin(x) - 0.25); });


    //matrix<double> A;
    //A = matrix<double>::randomDiagonal(4, -100, 100);
    //std::cout << "A:\n" << A << "\n";
    //matrix<double> T = (matrix<int>::random(4, 4, -100, 100));
    //A = T * A * (T.inverse_M());

    ////auto H = matrixfunction::sanitize_zeros(matrixfunction::hessenberg_upper_form(A), 1e-10);
    //std::cout << "A:\n" << A << "\n";
    ////std::cout << "T:\n" << T << "\n";
    ////std::cout << "hessenberg_form(A):\n" << H << "\n";
    ////std::cout << "A:\n" << A << "\n";

    //std::vector<std::complex<double>> Ans;
    //Ans=matrixfunction::compute_eigenvalues_3_qr(A);
    //for (auto i : Ans) {
    //    std::cout << i << "\n";
    //}
    


#endif


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
    {
        polynomial<double>a; a.outm_E = output_mode::FULL;
        polynomial<double>b;
        std::vector <double> roots{ 1,3,/*0.5,0.6,12 */};
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

        auto spl = Spline_interpolator<3,2, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << '\n' << spl << '\n';

        auto spl2 = Spline_interpolator<2, 1, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << '\n' << spl2 << '\n';

#define identifier (x-std::sin(x) - 0.25)
#define the_left_border -M_PI / 4
#define the_right_border M_PI / 4
        using SplineInterpolatorFunc = Spline<double>(*)(std::vector<std::pair<double, double>>);

        //x-std::sin(x) - 0.25
        counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<4, 3, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );  

        counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<2, 1, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<2, 1, double>",
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<1, 0, double>, // interpolator
            10, 20,                      // n, m_
            "Spline_interpolator<1,0, double>",
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        SHOW_INTERPOL_STAT(
            nuton_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            Lagrang_interpolation<double>,
            10, 10,
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            Alternativ_nuton_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            Alternativ_Lagrang_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );



    }
    
#elif AU_LOG == 6

std::vector<double> array = { 1,2,3,4,5,6 };
polynomial<double> polynom = w_k_T0(array, -1);

std::cout << "wk" << w_k_T0(array, -6).output_mode_set(output_mode::ABBREVIATED);
std::cout << "\nwx0" << (w_k_T0(array, 1))(1);

std::vector<std::pair<double, double>> array_xy = {
    {0.5, polynom(0.5) },{2.5, polynom(2.5)},
    {4.5, polynom(4.5) },{6.5, polynom(6.5)},
    {8.5, polynom(8.5) },{10.5, polynom(10.5)}
};


auto Func = [](double x) { return  1 + 10 * x + 10 * x * x + 10 * x * x * x + 10 * x * x * x * x; };
auto Array_xy = generatePoints_equally_sufficient_with_step_size(7, -4.0, 1.0, Func);

std::cout << "\nans:" << Lagrang_interpolation(Array_xy);
#elif AU_LOG == 7
    int Degree = 30;

    std::vector<double> array = { 1,2,3,4,5,6 };
    polynomial<double> polynom = w_k_T0(array, -1);



    std::cout << "wk" << w_k_T0(array, -6).output_mode_set(output_mode::ABBREVIATED);
    std::cout << "\nwx0" << (w_k_T0(array, 1))(1);

    std::vector<std::pair<double, double>> array_xy = {
        {0.5, polynom(0.5) },{2.5, polynom(2.5)},
        {4.5, polynom(4.5) },{6.5, polynom(6.5)},
        {8.5, polynom(8.5) },{10.5, polynom(10.5)}
    };

    array_xy=generatePoints_equally_sufficient_(Degree, -20., +20., polynom);
    Spline<double> a(-5,5,10);
    //std::cout << a;
    const uint64_t m=3, p=2;
    
    std::cout << "<m=" << m << ",p=" << p << ">=\n";
   

    draw_function<double>(Spline_interpolator<m, p>(array_xy));
    draw_function<double>(polynom);
    draw_function<double>(Spline_interpolator<m, p>(array_xy));


    
    //Spline_interpolator<3, 2>(array_xy);
    //std::cout << "<m=" << 4 << ",p=" << 1 << ">=\n";
    //Spline_interpolator<4, 1>(array_xy);
    //Spline_interpolator<4, 2>(array_xy);
    //Spline_interpolator<4, 3>(array_xy);
#elif AU_LOG == 8
#endif
    

   
    


    system("pause");
    return 0;
}
