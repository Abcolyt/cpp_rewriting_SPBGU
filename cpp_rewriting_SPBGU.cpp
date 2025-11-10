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
    counting_methods_2::Polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic(func,n,m+50, #func, __VA_ARGS__)


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
    std::vector<std::pair<T, T>> AddNoiseToPoints(
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
        auto noisy_points = AddNoiseToPoints(clean_points, 3, 0.2);
        noisy_points.insert(noisy_points.end(), clean_points.begin(), clean_points.end());
        std::sort(noisy_points.begin(), noisy_points.end());
        for (auto& I : noisy_points) {
            std::cout << "<" << I.first << ";" << I.second << ">\n";
        }
        std::cout << "\npolynomial function :\n" << "x*x + 20*x*x*x+ 444*x*x*x*x*x\n";

        std::cout << "\nNormal equations polynomial:\n" << counting_methods_2::aproximate::LeastSquaresNormal(noisy_points, degree_of_the_polynomial);
        std::cout << "\nOrthogonal polynomial:\n" << counting_methods_2::aproximate::LeastSquaresOrthogonal(noisy_points, degree_of_the_polynomial);
    }
    //

    void sem_4() {

        using namespace counting_methods_2::Polynomial_interpolation::nuton2;


#define PART_OF_TASK 4



#if PART_OF_TASK==0
        Z5_qr();
        // // // // 
        Z6_1(20, 7);

#define identifier (x-std::sin(x) - 0.25)
#define the_left_border  -M_PI / 4
#define the_right_border  M_PI / 4
        using SplineInterpolatorFunc = Spline<double>(*)(std::vector<std::pair<double, double>>);

        //x-std::sin(x) - 0.25
        counting_methods_2::Polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>,   // interpolator
            50, 50,                              // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );


        counting_methods_2::aproximate::output_of_characteristics_for_different_data_size_parameters::ShowAproximateStatistic([](double x) { return (x - std::sin(x) - 0.25); });

        // // // // // // // // //


#elif PART_OF_TASK == 4

        polynomial<double>a; a.outm_E = output_mode::FULL;
        polynomial<double>b;
        std::vector <double> roots{ 1,3,/*0.5,0.6,12 */ };
        a.the_root_constructor(roots);
        auto V = (polynomialfunctions::plnm_roots(a.get_first_derrivate(), DBL_EPSILON));
        std::cout << "polynomialfunctions:" << a << '\n';

        b = NutonInterpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\n:NutonInterpolation" << b;
        b = LagrangInterpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\nLagrang_interpolation:" << b;
        b = AlternativNutonInterpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
        b = polynomialfunctions::filter_large_epsilon(b, 1e-9);
        std::cout << "\nAlternativ_nuton_interpolation:" << b;
        b = AlternativLagrangInterpolation(generatePoints_equally_sufficient_(13, -12.0, 12.0, a));
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

        std::cout << "interpolinoms:NutonInterpolation, LagrangInterpolation, AlternativNutonInterpolation, AlternativLagrangInterpolation\n\n";


        SHOW_INTERPOL_STAT(
            NutonInterpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            LagrangInterpolation<double>,
            10, 20,
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            AlternativNutonInterpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );
        SHOW_INTERPOL_STAT(
            AlternativLagrangInterpolation<double>, // interpolator
            10, 20,                      // n, m_
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        std::cout << "spline: <3,2> <2,1> <1,0>\n";

        //x-std::sin(x) - 0.25
        counting_methods_2::Polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );

        counting_methods_2::Polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<2, 1, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<2, 1, double>",
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        counting_methods_2::Polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
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
        counting_methods_2::aproximate::output_of_characteristics_for_different_data_size_parameters::ShowAproximateStatistic([](double x) { return (x - std::sin(x) - 0.25); });
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
            abs_external_vector.push_back(std::abs(counting_methods_3::unsafe_integrate<TResult, TArg, integrate_method>(experimental_function, a, b, i, p_feature) - integral_of_the_function));
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
    namespace demonstration {
        namespace table {
            template<typename T>
            std::string anyToString(const T& value) {
                std::ostringstream oss;
                oss << value;
                return oss.str();
            }
            template<>std::string anyToString<std::string>(const std::string& value) {
                return value;
            }

            template<typename Tuple, typename Func, size_t... I>
            void applyWithIndexImpl(const Tuple& t, Func&& f, std::index_sequence<I...>) {
                (f(std::get<I>(t), I), ...);
            }

            template<typename Tuple, typename Func>
            void applyWithIndex(const Tuple& t, Func&& f) {
                applyWithIndexImpl(t, std::forward<Func>(f),
                    std::make_index_sequence<std::tuple_size_v<Tuple>>{});
            }

            template<typename... ColumnTypes>
            void printTableHeader(const std::vector<std::vector<std::string>>& headers,
                const std::vector<std::tuple<ColumnTypes...>>& data) {

                std::stringstream buffer;
                std::stringstream buffer_separator;

                constexpr size_t num_columns = sizeof...(ColumnTypes);

                std::vector<size_t> column_widths(num_columns, 0);

                for (size_t col = 0; col < num_columns; ++col) {
                    for (const auto& header_row : headers) {
                        if (col < header_row.size()) {
                            column_widths[col] = std::max(column_widths[col], header_row[col].size());
                        }
                    }
                }

                for (const auto& row : data) {
                    applyWithIndex(row, [&](const auto& element, size_t col_index) {
                        column_widths[col_index] = std::max(column_widths[col_index], anyToString(element).size());
                        });
                }
                for (auto& width : column_widths) {
                    width += 2;
                }

                buffer_separator << "+";
                for (size_t col = 0; col < num_columns; ++col) {
                    buffer_separator << std::string(column_widths[col], '-') << "+";
                }
                buffer_separator << "\n";
                buffer << buffer_separator.str();

                for (const auto& header_row : headers) {
                    buffer << "|";
                    for (size_t col = 0; col < num_columns; ++col) {
                        if (col < header_row.size()) {
                            buffer << std::setw(column_widths[col]) << std::left << header_row[col] << "|";
                        }
                        else {
                            buffer << std::string(column_widths[col], ' ') << "|";
                        }
                    }
                    buffer << "\n";
                }

                buffer << buffer_separator.str();

                for (const auto& row : data) {
                    buffer << "|";
                    applyWithIndex(row, [&](const auto& element, size_t col_index) {
                        buffer << std::setw(column_widths[col_index]) << std::left << anyToString(element) << "|";
                        });
                    buffer << "\n";
                }

                buffer << buffer_separator.str();
                std::cout << buffer.str();
            }
        }
        using namespace table;
        
        //Demonstration of the error of 4 methods:LEFT_RECTANGLE,MIDDLE_RECTANGLE,TRAPEZOID, SIMPSON
        template<typename TResult, typename TArg>
        void DispAbsBeatwen__4X(std::function<TResult(TArg)> experimental_function,
            TArg a, TArg b, TResult integral_of_the_function,
            counting_methods_3::PFeature p_feature = { 0,0 }, int n = 1) {
            std::vector<TResult> abs_external_vector1, abs_external_vector2, abs_external_vector3, abs_external_vector4;

            std::vector<std::vector<std::string>> headers = { {"(n)", "Error_LEFT_RECTANGLE","Error_MIDDLE_RECTANGLE","Error_TRAPEZOID", "Error_SIMPSON"} };
            std::vector<std::tuple<int, double, double, double, double>> data;
            for (int i = 1; i < n; i++)
            {
                try
                {
                    auto var1 = std::abs(counting_methods_3::unsafe_integrate<TResult, TArg, counting_methods_3::IntegrateMethod::LEFT_RECTANGLE>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var2 = std::abs(counting_methods_3::unsafe_integrate<TResult, TArg, counting_methods_3::IntegrateMethod::MIDDLE_RECTANGLE>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var3 = std::abs(counting_methods_3::unsafe_integrate<TResult, TArg, counting_methods_3::IntegrateMethod::NEWTON_COTES_2_POINT>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var4 = std::abs(counting_methods_3::unsafe_integrate<TResult, TArg, counting_methods_3::IntegrateMethod::SIMPSON>(experimental_function, a, b, i, p_feature) - integral_of_the_function);

                    data.push_back(std::tuple{ i, var1,var2,var3,var4 });

                    abs_external_vector1.push_back(var1);
                    abs_external_vector2.push_back(var2);
                    abs_external_vector3.push_back(var3);
                    abs_external_vector4.push_back(var4);
                    std::cout << "continued i:" << i << "\n";
                }
                catch (const std::exception&)
                {
                    continue;
                }

            }
            std::vector<std::function<double(double)>> funcs;

            auto f_abs_n1 = [&abs_external_vector1](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector1.size()) {
                    return abs_external_vector1[index];
                }
                return 0.0;
                };
            funcs.push_back(f_abs_n1);
            auto f_abs_n2 = [&abs_external_vector2](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector2.size()) {
                    return abs_external_vector2[index];
                }
                return 0.0;
                };
            funcs.push_back(f_abs_n2);
            auto f_abs_n3 = [&abs_external_vector3](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector3.size()) {
                    return abs_external_vector3[index];
                }
                return 0.0;
                };
            funcs.push_back(f_abs_n3);
            auto f_abs_n4 = [&abs_external_vector4](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector4.size()) {
                    return abs_external_vector4[index];
                }
                return 0.0;
                };
            funcs.push_back(f_abs_n4);

            sem_5::demonstration::printTableHeader(headers, data);

#if __has_include(<SFML/Graphics.hpp>)
            draw_functions(funcs, 0.0, static_cast<double>(n + 1));

#else
            std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
                << "not working draw_functions(functions_normal)\n"
                << "not working draw_functions(functions_ortog)\n";
#endif
        }

    }
    using namespace demonstration;
    
    namespace integral {

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

    TestConfig<double, double> exemple_1{
           [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
           [](double x) { return (-109680 * sin((9 * x) / 2) - 3948480 * cos((9 * x) / 2) + 1062243 * exp((19 * x) / (24)) * sin((7 * x) / 2) + 202332 * exp((19 * x) / (24)) * cos((7 * x) / 2)) / (2963645 * exp(x / 8)) + (5 * x * x) / 2 + 3746148 / (2963645); },
           0.7, 3.2,
           20.235344669492576,
           1e-9,
           30,
           "variant_10_alpha_beta_0",
           0.0,0.0,
           1
    };
    TestConfig<double, double> example_2{
           [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
           [](double x) { return (-109680 * sin((9 * x) / 2) - 3948480 * cos((9 * x) / 2) + 1062243 * exp((19 * x) / (24)) * sin((7 * x) / 2) + 202332 * exp((19 * x) / (24)) * cos((7 * x) / 2)) / (2963645 * exp(x / 8)) + (5 * x * x) / 2 + 3746148 / (2963645); },
           7.0/10, 3.2,
           24.142092678549481,
           1e-9,
           30,
           "variant_10_alpha_beta_0",
           0.0,.25,
           1
    };
    
    }

    void sem_5_part1() {
        using namespace integral;
        auto configurate = exemple_1;
#define EXECUTION_PART 10

#if EXECUTION_PART == 10 || EXECUTION_PART == all
        DispAbsBeatwen__4X<double, double>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 800);
#endif
        ////Part No. 1: Quadrature formulas of Newton-Kot(e)ca and Gauss
        //1.1 calculations of a certain integral using compound quadrature formulas 
        configurate = example_2;
#if EXECUTION_PART == 11 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::LEFT_RECTANGLE>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::MIDDLE_RECTANGLE>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::TRAPEZOID>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::SIMPSON>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
#endif

        //Methods for calculating a definite integral using compound quadrature formulas
        //based on 3 - point Newton - Kot formulas(e)ca and Gauss.
#if EXECUTION_PART == 12 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_3_POINT> (configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_4_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::NEWTON_COTES_5_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
#endif
        //1.3. Graph of the dependence of the absolute error on the number of splits 
        //of the integration interval of each quadrature formula from paragraphs 1.1 - 2.
#if EXECUTION_PART == 13 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::GAUS_3_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, counting_methods_3::IntegrateMethod::GAUS_4_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 20);
#endif
        ////Part No. 2: Methods for estimating the error of composite quadrature formulas
        //A definite integral with a given accuracy 𝜀 using the composite 3 - point quadrature formula of Newton - Cot(e)ca.
        // 
        //Estimation of the error by the Richardson method.
        //The rate of convergence according to Aitken's rule.
        // 
        //The length of the step ℎ of the partition.
        counting_methods_3::IntegrationResult<double> result;
#if EXECUTION_PART == 21 || EXECUTION_PART == all
        result = counting_methods_3::adaptive_integrate<double, double,
            counting_methods_3::IntegrateMethod::NEWTON_COTES_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, 1, { configurate.alpha, configurate.beta });
        std::cout <<"IntegrateMethod::NEWTON_COTES_3_POINT\n" << result;
#endif

        // N2.1 but using 3-point Gauss formulas:
#if EXECUTION_PART == 22 || EXECUTION_PART == all
        result = counting_methods_3::adaptive_integrate<double, double,
            counting_methods_3::IntegrateMethod::GAUS_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, 1, { configurate.alpha, configurate.beta });
        std::cout <<"IntegrateMethod::GAUS_3_POINT\n" << result;
#endif

        //Using the Aitken convergence rate estimate, select the h_opt step for either formulas 2.1 or 2.2.
        // 
        //Сompare it with the step calculated in 2.1 or 2.2.
#if EXECUTION_PART == 23 || EXECUTION_PART == all
        auto simple_result = counting_methods_3::optimized_adaptive_integration<double, double,
            counting_methods_3::IntegrateMethod::GAUS_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, { configurate.alpha, configurate.beta });
        std::cout << "IntegrateMethod::GAUS_3_POINT\n" << simple_result;
#endif
        ////
#undef EXECUTION_PART
    }
    namespace differential_equation {



    template<typename Targ, typename Tresult>
    struct TestConfig {
        Targ c2;
        Targ x0;
        Tresult y0;
        Targ x_max;
        const Targ h;

        std::function<Tresult(Targ, Tresult)> f;
        std::function<Tresult(Targ)> exact_y;

        double global_accuracy;
        double local_accuracy;
    };

    TestConfig<double, double> basic{
        0.5,
        0.0,
        1.0,
        1.0,
        0.1,
        [](double x, double y) {return x * y; },
        [](double x) {return std::exp(x * x / 2); },
        1e-6,
        1e-6
    };
    
    }


    template<typename Targ, typename Tresult>
    std::vector<std::pair<Targ, Tresult>> runge_kutta_2nd_order(
        std::function<Tresult(Targ, Tresult)> f,  // Right-hand side function f(x, y)
        /*Tresult(*exact_y)(Targ x),*/   // Exact solution y(x) for error analysis
        const Targ x0,                 // Initial x
        const Tresult y0,                // Initial y
        Targ h,                  // Step size
        Targ x_max,              // Maximum x value
        Targ c2                  // Parameter ξ for RK2 scheme
    ) {
        // coefficients
        Targ a21 = c2;
        Targ b2 = static_cast<Targ>(1) / (2 * c2);
        Targ b1 = 1 - static_cast<Targ>(1) / (2 * c2);

        //Ans
        std::vector<std::pair<Targ, Tresult>> solution;
        solution.push_back({ x0, y0 });


        Targ x_current = x0;
        Tresult y_current = y0;

        while (x_current < x_max) {
            // RK2 stages
            Tresult k1 = h * f(x_current, y_current);
            Tresult k2 = h * f(x_current + c2 * h, y_current + a21 * k1);

            // Update solution
            y_current = y_current + b1 * k1 + b2 * k2;
            x_current += h;

            solution.push_back({ x_current, y_current });
        }

        return solution;
    }
    template<typename Targ, typename Tresult>
    std::vector<std::pair<Targ, Tresult>> runge_kutta_2nd_order_wrapper(differential_equation::TestConfig<Targ,Tresult> conf) {
        return runge_kutta_2nd_order<Targ, Tresult>(conf.f, conf.x0, conf.y0, conf.h, conf.x_max, conf.c2);
    }
    void sem_5_part2(){
#define EXECUTION_PART 11
        /*
        The Cauchy problem :
            dy1(𝑥) / dx = 𝐴_𝑦2(𝑥),
            dy2(𝑥) / dx = -B_𝑦1(𝑥),
            y1(0) = Bπ,
            y2(0) = Aπ,
        */

        
        ////Part No. 1: Calculation schemes of the Runge-Kutta method with a constant step
        //A second-order calculation scheme based on second-order conditions for the 2-stage explicit Rungi-Kutta method for the c2 parameter
#if EXECUTION_PART == 11 || EXECUTION_PART == all
         

        // Function definitions (replace with your actual functions)
        //auto f = [](T x,T2 y){ return x * y; }; // Right-hand side: dy/dx = f(x, y)
        //auto exact_y = [](T x){ return exp(x * x / 2); }; // Exact solution for error analysis



#if 0
        auto runge_abs_error = [&y1,&y2](s) {return (y1-y2)/(std::pow(2,s) -1); };
#endif
        // Implementation of RK scheme construction with parameter c2 = ξ
        auto solution = runge_kutta_2nd_order_wrapper<double, double>(differential_equation::basic);

        for(auto& var :solution)
        {
            std::cout << "{" << var.first << ";" << var.second << "} exact sol:{"<< var.first<<";"<< differential_equation::basic.exact_y(var.first) << "}"
                <<"difference:"<< var.second - differential_equation::basic.exact_y(var.first)<<"\n";
        }

#endif  

        // 1.2 Implement constant-step RK method with total error estimation using Runge method
        //     (ε = 1e-4), initial step selection according to algorithm
#if EXECUTION_PART == 12 || EXECUTION_PART == all
#endif

        ////

        ////Part No. 2: Calculation schemes of the Runge-Kutta method with a constant step
        // 2.1 Implement automatic step-size control using 2-stage RK 2nd order method
        //     with local error estimation (ρ = 1e-5) and Runge method
#if EXECUTION_PART == 21 || EXECUTION_PART == all

#endif
        ////

        ////Part No. 3: Calculation schemes of the Runge-Kutta method with a constant step
        // 3.1 Implement constant-step and automatic-step methods using classical RK schemes
        //     of 3rd or 4th order (opponent scheme)
#if EXECUTION_PART == 31 || EXECUTION_PART == all
#endif

        // 3.2 Determine integration step h for constant-step methods (2-stage RK 2nd order
        //     and opponent scheme) that provides solution with accuracy ε = 1e-4
        //     Plot true total error vs x
#if EXECUTION_PART == 32 || EXECUTION_PART == all
#endif

        // 3.3 For automatic step-size control methods (2-stage RK 2nd order and opponent):
        // Analyze reliability and efficiency of implemented algorithms
        // 3.3.1 Plot step size vs x
#if EXECUTION_PART == 331 || EXECUTION_PART == all
#endif

        // 3.3.2 Plot ratio of true local error to estimated local error vs x
#if EXECUTION_PART == 332 || EXECUTION_PART == all
#endif

        // 3.3.3 Plot number of right-hand side evaluations vs accuracy ε
#if EXECUTION_PART == 333 || EXECUTION_PART == all
#endif
        ////

#undef EXECUTION_PART

    }
}

/////////////////////////////////////////////////////////


int main() {

    

    sem_4::sem_4();
    sem_5::sem_5_part1();
    //sem_5::sem_5_part2();



    system("pause");
    return 0;
}



