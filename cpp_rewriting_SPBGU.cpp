#include "file_h\project_libraries.h"


//uncomment it and you can use it carefully. 
//Performance is not guaranteed
namespace sem_3 {
    void sem_3() {

#if 1
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
    }

}
//Assignments for 4 semesters were working at the time of completion. 
//It should work in the vast majority of cases.
//Topics: interpolation, eigenvalues
namespace sem_4 {

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
        auto T = (lup.P * b); T.set_output_mode(output_mode::ABBREVIATED);
        std::cout << "b:\n" << T.set_output_mode(output_mode::ABBREVIATED) << "\n";
        matrix<double> LU = (lup.L * lup.U); LU.set_output_mode(output_mode::ABBREVIATED);
        /*std::cout << "l:\n" << ((LU.get_output_mode()) == output_mode::FULL) << "\n";*/
        matrix<double> PA = (lup.P * A); PA.set_output_mode(output_mode::ABBREVIATED);
        std::cout << "L * U:\n" << LU << "\nP * A:\n" << PA << std::endl;

    }

    void Z5_2(int Size = 4) {
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
        matrix<double> A = matrix<double>::zeros(4, 4);
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
        for (auto& i : A_.determinant().plnm_roots())
        {
            j++;
            std::cout << "root_[" << j << "]= " << i.first << "\n";
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
            std::cout << "eig[" << j << "]= " << i << "\n";
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
    void Z6_1(int size, int degree_of_the_polynomial) {
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
        std::cout << "\npolynomial function :\n" << "x*x + 20*x*x*x+ 444*x*x*x*x*x\n";

        std::cout << "\nNormal equations polynomial:\n" << counting_methods_2::aproximate::least_squares_normal(noisy_points, degree_of_the_polynomial);
        std::cout << "\nOrthogonal polynomial:\n" << counting_methods_2::aproximate::least_squares_orthogonal(noisy_points, degree_of_the_polynomial);
    }
    //

    void sem_4() {

        using namespace counting_methods_2::Polynomial_interpolation::nuton2;


#define PART_OF_TASK 0



#if PART_OF_TASK==0
        Z5_qr();
        // // // // 
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


#elif PART_OF_TASK == 4

        polynomial<double>a; a.outm_E = output_mode::FULL;
        polynomial<double>b;
        std::vector <double> roots{ 1,3,/*0.5,0.6,12 */ };
        a.the_root_constructor(roots);
        auto V = (polynomialfunctions::plnm_roots(a.get_first_derrivate(), DBL_EPSILON));
        std::cout << "polynomialfunctions:" << a << '\n';

        b = nuton_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\n:nuton_interpolation" << b;
        b = Lagrang_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\nLagrang_interpolation:" << b;
        b = Alternativ_nuton_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\nAlternativ_nuton_interpolation:" << b;
        b = Alternativ_Lagrang_interpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\nAlternativ_Lagrang_interpolation:" << b;

        auto spl4 = Spline_interpolator<4, 3, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << "\nSpline_interpolator<4,3, double>:" << spl4 << "\n";

        auto spl3 = Spline_interpolator<3, 2, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << "\nSpline_interpolator<3,2, double>:" << spl3 << "\n";

        auto spl2 = Spline_interpolator<2, 1, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << "\nSpline_interpolator<2, 1, double>:" << spl2 << "\n";


        auto spl1 = Spline_interpolator<1, 0, double>(generatePoints_equally_sufficient_(5, -12.0, 12.0, a));
        std::cout << "\nSpline_interpolator<1, 0, double>:" << spl1 << "\n";
#define identifier (x-std::sin(x) - 0.25)
#define the_left_border -M_PI / 4
#define the_right_border M_PI / 4
        using SplineInterpolatorFunc = Spline<double>(*)(std::vector<std::pair<double, double>>);

        std::cout << "interpolinoms:nuton_interpolation, Lagrang_interpolation, Alternativ_nuton_interpolation, Alternativ_Lagrang_interpolation\n\n";


        SHOW_INTERPOL_STAT(
            nuton_interpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            Lagrang_interpolation<double>,
            10, 20,
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

        std::cout << "spline: <3,2> <2,1> <1,0>\n";

        //x-std::sin(x) - 0.25
        counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>, // interpolator
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



#elif PART_OF_TASK==5
        Z5_2(4);
        Z5_1_2_const();
        Z5_plus();
        Z5_qr();

#elif PART_OF_TASK == 6
        Z6_1(20, 7);
        counting_methods_2::aproximate::show_aproximate_statistic([](double x) { return (x - std::sin(x) - 0.25); });
#endif




    }
}

namespace sem_5 {
    template<typename TResult, typename TArg, counting_methods_3::IntegrateMethod integrate_method>
    void DispAbsBeatwen(std::function<TResult(TArg)> experimental_function,
        TArg a, TArg b, TResult integral_of_the_function,
        counting_methods_3::PFeature p_feature={0,0},int n=1) {
        std::vector<TResult> abs_external_vector;
        for (int i = 0; i < n; i++)
        {
            try
            {
            abs_external_vector.push_back(std::abs(counting_methods_3::integrate<TResult, TArg, integrate_method>(experimental_function, a, b, i, p_feature) - integral_of_the_function));
            std::cout << "continued i:" << i<<"\n";
            }
            catch (const std::exception&)
            {
                continue;
            }

        }
        auto f_abs_n = [&abs_external_vector](double x) -> double {
            int index = static_cast<int>(std::floor(x)); 
            if (index >= 0 && index < abs_external_vector.size()) {
                return abs_external_vector[index];
            }
            return 0.0; 
            };
        std::vector<std::function<double(double)>> funcs;
        funcs.push_back(f_abs_n);
#if __has_include(<SFML/Graphics.hpp>)
        draw_functions(funcs, 0.0, static_cast<double>(n + 1));

#else
        std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
            << "not working draw_functions(functions_normal)\n"
            << "not working draw_functions(functions_ortog)\n";
#endif
        //std::cout << "end ALL\n";
        //draw_functions(std::vector{ f_abs_n }, 0,0, static_cast<double>(n + 1));
    }

    //show all the absolute values of the differences with the reference
    template<typename TResult, typename TArg>
    void DispAllAbs(std::function<TResult(TArg)> experimental_function, TArg a, TArg b, TResult integral_of_the_function, counting_methods_3::PFeature p_feature = { 0,0 }, int n = 1) {
        for (int i=0; i< static_cast<int>(counting_methods_3::IntegrateMethod::LEFT_RECTANGLE); i++)
        {
            DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(i)>(experimental_function, a,b,integral_of_the_function, p_feature, n);
        }
    }
    template<typename TArg, typename TResult>
    struct TestConfig {
        std::function<TResult(TArg)> function;      // Интегрируемая функция
        std::function<TResult(TArg)> reference;     // Первообразная (эталонная функция)
        TArg a;
        TArg b;
        TResult expected;                           // Ожидаемый результат
        TResult tolerance;                          // Допустимая погрешность
        int points;
        std::string name;                           // Имя конфигурации( для отладки)
        TArg alpha;                                 // Дополнительный параметр
        TArg beta;
        int max_iteration = 10;
    };

    TestConfig<double, double> conf{
           [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
           [](double x) { return (-109680 * sin((9 * x) / 2) - 3948480 * cos((9 * x) / 2) + 1062243 * exp((19 * x) / (24)) * sin((7 * x) / 2) + 202332 * exp((19 * x) / (24)) * cos((7 * x) / 2)) / (2963645 * exp(x / 8)) + (5 * x * x) / 2 + 3746148 / (2963645); },
           0.7, 3.2,
           20.235345,
           1e-9,
           30,
           "variant_10_alpha_beta_0",
           0.0,0.0,
           1
    };
    TestConfig<double, double> conf2{
           [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
           [](double x) { return (-109680 * sin((9 * x) / 2) - 3948480 * cos((9 * x) / 2) + 1062243 * exp((19 * x) / (24)) * sin((7 * x) / 2) + 202332 * exp((19 * x) / (24)) * cos((7 * x) / 2)) / (2963645 * exp(x / 8)) + (5 * x * x) / 2 + 3746148 / (2963645); },
           0.7, 3.2,
           24.142092678433,
           1e-9,
           30,
           "variant_10_alpha_beta_0",
           0.0,1.0/4,
           1
    };
    TestConfig<double, double> conf3{
               [](double x) { return 1.5 * cos(3.7 * x) * exp(4 * x / 7) + 3 * sin(2.5 * x) * exp(3 * x / 4) + 3 * x;  },
               [](double x) { return (27195 * exp((4 * x) / 7) * sin((3.7 * x))) / (68681) + (4200 * exp((4 * x) / 7) * cos(3.7 * x)) / (68681) + (36 * exp(0.75 * x) * sin(2.5 * x)) / (109) - (120 * exp(0.75 * x) * cos(2.5 * x)) / (109) + (3 * x * x) / 2 + 7783920 / (7486229); },
               1.5, 3.0,
               5.6085702,
               1e-9,
               30,
               "variant_20_alpha_0_beta_5__6",
               0,5.0 / 6,
               3
    };
    void sem_5() {
        counting_methods_3::PFeature p={ 0,0};
        auto f = [](double x) {return 13*std::cos(3.5*x) *std::exp(2*x/ 3) + 6*std:: sin(4.5*x) *exp(x /8) + 5*x; };
        
        //((e ^ 2) ^ (1 / 5) * (103296 * sin(144) - 37186560 * cos(144)) + (e ^ 2) ^ (1 / 15) * (103554360 * e ^ 2 * sin(112) + 1972464 * e ^ 2 * cos(112)) + (e ^ 7) ^ (1 / 80) * (37186560 * cos(63 / 2) - 103296 * sin(63 / 2)) + (e ^ 7) ^ (1 / 15) * (-103554360 * sin(49 / 2) - 1972464 * cos(49 / 2)) + 6798220455) / (278901352) = 21.92472

        auto configurate_alp_bet_eq0 = conf,configurate = conf2;



        auto result = counting_methods_3::adaptive_integrate_richardson<double, double,
            counting_methods_3::IntegrateMethod::NEWTON_COTES_3_POINT>(
                configurate.function, configurate.a, configurate.b, 1e-6, 1, { 0, 0 });

        std::cout << "RESULT STATISTIC:" << std::endl;
        std::cout << "integral value : " << result.value << std::endl;
        std::cout << "refined_value: " << result.refined_value << std::endl;
        std::cout << "estimated_error: " << result.estimated_error << std::endl;
        std::cout << "step_size: " << result.step_size << std::endl;
        std::cout << "final_intervals number: " << result.final_intervals << std::endl;

        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::GAUS_3_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::GAUS_4_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 20);


        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(0)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(1)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(2)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(3)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(4)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(5)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);
        DispAbsBeatwen<double, double, static_cast<counting_methods_3::IntegrateMethod>(6)>(configurate_alp_bet_eq0.function, configurate_alp_bet_eq0.a, configurate_alp_bet_eq0.b, 20.235345, { configurate_alp_bet_eq0.alpha,configurate_alp_bet_eq0.beta }, 100);




        //counting_methods_3::integrate<double,double,counting_methods_3::IntegrateMethod::LEFT_RECTANGLE>(f,p.a,p.b);


        //const std::vector<float> execution_logic={1.1,1.2,1.3};

        ////Part No. 1: Quadrature formulas of Newton-Kot(e)ca and Gauss
        //1.1 calculations of a certain integral using compound quadrature formulas 
        // Beginning:

        //ens;

        //Methods for calculating a definite integral using compound quadrature formulas
        //based on 3 - point Newton - Kot formulas(e)ca and Gauss.
        // Beginning:

        //ens;

        //1.3. Graph of the dependence of the absolute error on the number of splits 
        //of the integration interval of each quadrature formula from paragraphs 1.1 - 2.
        // Beginning:

        //ens;

        ////

        ////Part No. 2: Methods for estimating the error of composite quadrature formulas
        //A definite integral with a given accuracy 𝜀 using the composite 3 - point quadrature formula of Newton - Cot(e)ca.
        // 
        //Estimation of the error by the Richardson method.
        //The rate of convergence according to Aitken's rule.
        // 
        //The length of the step ℎ of the partition.
        // Beginning:

        //ens;

        // N2.1 but using 3-point Gauss formulas:
        // Beginning:

        //ens;

        //Using the Aitken convergence rate estimate, select the h_opt step for either formulas 2.1 or 2.2.
        // 
        //Сompare it with the step calculated in 2.1 or 2.2.
        // Beginning:

        //ens;
        ////
    }
}
#if 0
template<typename IntegratedFunction>
inline double ApproximateValueOfTheIntegral(matrix<double> nodes_in_the_integration_gap, IntegratedFunction F, double a, double b) {
    if (interval_partitioning_scheme.is_vector() == false || (left_boundary_of_the_integration_gap > right_boundary_of_the_integration_gap) || )


    {
        throw std::exception("interval_partitioning_scheme is not a vector ");
    }
    if (interval_partitioning_scheme.is_vertical_vector()) {
        interval_partitioning_scheme = interval_partitioning_scheme.transpose();
    }
    double ans = 0;
    matrix<double> coeff = counting_methods_3::СoefficNewtonCotes(nodes_in_the_integration_gap, a, b);

    for (int i = 0; i < 4; i++)
    {
        ans += F(nodes_in_the_integration_gap[i][]);
    }

    return ans * (b - a);
}

#endif
/////////////////////////////////////////////////////////


#if 0
template<typename P>
std::vector<std::pair<std::complex<P>, int>> analytical_polynomial_roots(const polynomial<P>& plnm) {
    std::vector<std::pair<std::complex<P>, int>> roots;
    size_t deg = plnm.get_deg();

    if (deg > 5) {
        return roots; 
    }

    std::vector<std::complex<P>> complex_roots;

    switch (deg-1) {
    case 0:
        break;
    case 1:
        complex_roots.push_back(std::complex<P>(-plnm[0] / plnm[1], 0));
        break;
    case 2:
        complex_roots = solve_quadratic(plnm[2], plnm[1], plnm[0]);
        break;
    case 3:
        complex_roots = solve_cubic(plnm[3], plnm[2], plnm[1], plnm[0]);
        break;
    case 4:
        complex_roots = solve_quartic(plnm[4], plnm[3], plnm[2], plnm[1], plnm[0]);
        break;
    }

    for (const auto& root : complex_roots) {
        if (std::abs(root.imag()) < 1e-10) {
            roots.push_back(std::make_pair(root.real(), 0));
        }
        else
            roots.push_back(std::make_pair(root,0));
    }


    return roots;
}
/////////////
// Тестовая функция
template<typename P>
void test_analytical_roots() {
    std::cout << "=== Testing Analytical Polynomial Roots ===\n\n";

    // Тест 1: Полином 1-й степени: 2x - 4 = 0 (корень: 2)
    {
        polynomial<P> p1;
        p1.newsize(2);
        p1[0] = -4;
        p1[1] = 2;

        std::cout << "1st degree polynomial: 2x - 4 = 0\n";
        auto roots = analytical_polynomial_roots(p1);
        for (const auto& root : roots) {
            std::cout << "Root: " << root.first << ", Iterations: " << root.second << std::endl;
        }
        std::cout << "Expected: Root = 2.0, Iterations = 0\n\n";
    }

    // Тест 2: Полином 2-й степени: x² - 5x + 6 = 0 (корни: 2, 3)
    {
        polynomial<P> p2;
        p2.newsize(3);
        p2[0] = 6;
        p2[1] = -5;
        p2[2] = 1;

        std::cout << "2nd degree polynomial: x^2 - 5x + 6 = 0\n";
        auto roots = analytical_polynomial_roots(p2);
        for (const auto& root : roots) {
            std::cout << "Root: " << root.first << ", Iterations: " << root.second << std::endl;
        }
        std::cout << "Expected: Roots =~= 2.0, 3.0, Iterations = 0\n\n";
    }

    // Тест 3: Полином 3-й степени: x³ - 6x² + 11x - 6 = 0 (корни: 1, 2, 3)
    {
        polynomial<P> p3;
        p3.newsize(4);
        p3[0] = -6;
        p3[1] = 11;
        p3[2] = -6;
        p3[3] = 1;

        std::cout << "3rd degree polynomial: x^3 - 6x^2 + 11x - 6 = 0\n";
        auto roots = analytical_polynomial_roots(p3);
        for (const auto& root : roots) {
            std::cout << "Root: " << root.first << ", Iterations: " << root.second << std::endl;
        }
        std::cout << "Expected: Roots =~= 1.0, 2.0, 3.0, Iterations = 0\n\n";
    }

    // Тест 4: Полином 4-й степени: x⁴ - 10x³ + 35x² - 50x + 24 = 0 (корни: 1, 2, 3, 4)
    {
        polynomial<P> p4;
        p4.newsize(5);
        p4[0] = 24;
        p4[1] = -50;
        p4[2] = 35;
        p4[3] = -10;
        p4[4] = 1;

        std::cout << "4th degree polynomial: x^4 - 10x^3 + 35x^2 - 50x + 24 = 0\n";
        auto roots = analytical_polynomial_roots(p4);
        for (const auto& root : roots) {
            std::cout << "Root: " << root.first << ", Iterations: " << root.second << std::endl;
        }
        std::cout << "Expected: Roots =~= 1.0, 2.0, 3.0, 4.0, Iterations = 0\n\n";
    }

    // Тест 5: Полином 2-й степени с комплексными корнями: x² + 1 = 0 (вещественных корней нет)
    {
        polynomial<P> p5;
        p5.newsize(3);
        p5[0] = 1;
        p5[1] = 0;
        p5[2] = 1;

        std::cout << "2nd degree polynomial with complex roots: x^2 + 1 = 0\n";
        auto roots = analytical_polynomial_roots(p5);
        std::cout << "Found " << roots.size() << " real roots\n";
        for (const auto& root : roots) {
            std::cout << "Root: " << root.first << ", Iterations: " << root.second << std::endl;
        }
        std::cout << "Expected: No real roots\n\n";
    }

    // Тест 6: Полином 5-й степени (должен вернуть пустой вектор)
    {
        polynomial<P> p6;
        p6.newsize(6);
        p6[0] = 1; p6[1] = 2; p6[2] = 3; p6[3] = 4; p6[4] = 5; p6[5] = 6;

        std::cout << "5th degree polynomial (should use numerical method)\n";
        auto roots = analytical_polynomial_roots(p6);
        std::cout << "Found " << roots.size() << " roots (analytical method not applicable)\n";
        std::cout << "Expected: Empty vector (analytical method only for degrees < 5)\n\n";
    }
}
#endif

int main() {
#if 1
   
    sem_5::sem_5();

#else
    using namespace counting_methods_3;
    using namespace sem_5;  

    //std::cout << special::ComputeBetaZeroIntegral<double>(1., 10., 1., 1., 2)<<"\n";
    // 
    //std::cout << "First"  << std::fixed << std::setprecision(12) << special::integrated_function<double>(1., 7., 0., 8., 0., 0.99, 2) << "\n";
    //std::cout << "Second" << std::fixed << std::setprecision(12) << special::integrated_function<double>(1., 7., 1., 0., 0.7, 0., 2) << "\n";
    //First52.888989776304
    //Second48.299771618350

    std::cout << "Second" /*<< std::fixed << std::setprecision(12)*/ <<special::incomplete_beta(0.3,2.,3.) << "\n";
    auto configurate = conf2;
    



    std::cout << "Second" /*<< std::fixed << std::setprecision(30) */<< special::IntegralManager<double>(1., 7., 1., 8., 0.7, 0.99, 2)<<"\n";
    
    std::cout << "NEWTON_COTES_3_POINT, ANS:" << integrate<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_3_POINT>(configurate.function, configurate.a, configurate.b, 10'000, { configurate.alpha,configurate.beta}) << "end\n";

    std::cout << "GAUS_3_POINT, ANS:" << integrate<double, double, counting_methods_3::IntegrateMethod::GAUS_3_POINT>(configurate.function, configurate.a, configurate.b, 100, { configurate.alpha,configurate.beta }) << "end\n";

    //sem_5::sem_5();
    //GetNGausCoefficientWithP_x_FunctionConstants(3);
    //GetNGausCoefficientWithP_x_FunctionConstants(3, 0.7, 3.2, 0, 0);
    //auto Tay = conf2;
    //std::cout << "ANS:" << integrate<double, double, counting_methods_3::IntegrateMethod::GAUS_3_POINT>(Tay.function, Tay.a, Tay.b, 30, { Tay.alpha,Tay.beta }) << "end\n";

   // std::cout << "ANS:" << integrate<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_5_POINT>(conf.function, conf.a, conf.b, conf.points, { conf.alpha,conf.beta })<<"end\n";
    //std::cout << "ANS:" << integrate<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_3_POINT>(conf2.function, conf2.a, conf2.b, 1'00'000'000, { conf2.alpha,conf2.beta }) << "end\n";

#endif



    system("pause");
    return 0;
}



