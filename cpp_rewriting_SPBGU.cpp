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
    counting_methods_2::polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic(func,n,m+50, #func, __VA_ARGS__)


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


        std::vector<std::pair<double, double>> clean_points = counting_methods_2::polynomial_interpolation::nuton2::generatePoints_equally_sufficient_(size, the_left_border, the_right_border, [](double x) { return identifier; });
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

        using namespace counting_methods_2::polynomial_interpolation::nuton2;


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
        counting_methods_2::polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
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
        counting_methods_2::polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );

        counting_methods_2::polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<2, 1, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<2, 1, double>",
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        counting_methods_2::polynomial_interpolation::nuton2::output_of_characteristics_for_different_data_size_parameters::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
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
    template<typename TResult, typename TArg, IntegrateMethod integrate_method>
    void DispAbsBeatwen(std::function<TResult(TArg)> experimental_function,
        TArg a, TArg b, TResult integral_of_the_function,
        PFeature p_feature={0,0},int n=1) {
        std::vector<TResult> abs_external_vector;
        for (int i = 0; i < n; i++)
        {
            try
            {
            abs_external_vector.push_back(std::abs(unsafe_integrate<integrate_method, TResult, TArg>(experimental_function, a, b, i, p_feature) - integral_of_the_function));
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
        DrawFunctions(funcs, 0.0, static_cast<double>(n + 1));

#else
        std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
            << "not working DrawFunctions(functions_normal)\n"
            << "not working DrawFunctions(functions_ortog)\n";
#endif
        //std::cout << "end ALL\n";
        //DrawFunctions(std::vector{ f_abs_n }, 0,0, static_cast<double>(n + 1));
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
            PFeature p_feature = { 0,0 }, int n = 1) {
            std::vector<TResult> abs_external_vector1, abs_external_vector2, abs_external_vector3, abs_external_vector4;

            std::vector<std::vector<std::string>> headers = { {"(n)", "Error_LEFT_RECTANGLE","Error_MIDDLE_RECTANGLE","Error_TRAPEZOID", "Error_SIMPSON"} };
            std::vector<std::tuple<int, double, double, double, double>> data;
            for (int i = 1; i < n; i++)
            {
                try
                {
                    auto var1 =std::abs(unsafe_integrate<IntegrateMethod::LEFT_RECTANGLE, TResult, TArg>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var2 = std::abs(unsafe_integrate<IntegrateMethod::MIDDLE_RECTANGLE, TResult, TArg>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var3 = std::abs(unsafe_integrate<IntegrateMethod::NEWTON_COTES_2_POINT, TResult, TArg>(experimental_function, a, b, i, p_feature) - integral_of_the_function),
                        var4 = std::abs(unsafe_integrate<IntegrateMethod::SIMPSON, TResult, TArg>(experimental_function, a, b, i, p_feature) - integral_of_the_function);

                    data.push_back(std::tuple{ i, var1,var2,var3,var4 });

                    abs_external_vector1.push_back(var1);
                    abs_external_vector2.push_back(var2);
                    abs_external_vector3.push_back(var3);
                    abs_external_vector4.push_back(var4);
                    //std::cout << "continued i:" << i << "\n";
                }
                catch (const std::exception&)
                {
                    continue;
                }

            }
            std::vector<std::function<double(double)>> funcs;
#if 1
            auto f_abs_n1 = createLambda(abs_external_vector1);
#else
            auto f_abs_n1 = [&abs_external_vector1](double x) -> double {
                            int index = static_cast<int>(std::floor(x));
                            if (index >= 0 && index < abs_external_vector1.size()) {
                                return abs_external_vector1[index];
                            }
                            return 0.0;
                            };
#endif    
            funcs.push_back(f_abs_n1);

            auto f_abs_n2 = createLambda(abs_external_vector2);
            
            /*[&abs_external_vector2](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector2.size()) {
                    return abs_external_vector2[index];
                }
                return 0.0;
                };*/
            funcs.push_back(f_abs_n2);

            auto f_abs_n3 = createLambda(abs_external_vector3);
            /*[&abs_external_vector3](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector3.size()) {
                    return abs_external_vector3[index];
                }
                return 0.0;
                };*/
            funcs.push_back(f_abs_n3);
            auto f_abs_n4 = createLambda(abs_external_vector4);
            
            /*[&abs_external_vector4](double x) -> double {
                int index = static_cast<int>(std::floor(x));
                if (index >= 0 && index < abs_external_vector4.size()) {
                    return abs_external_vector4[index];
                }
                return 0.0;
                };*/
            funcs.push_back(f_abs_n4);

            sem_5::demonstration::printTableHeader(headers, data);

#if __has_include(<SFML/Graphics.hpp>)
            DrawFunctions(funcs, 0.0, static_cast<double>(n + 1));

#else
            std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
                << "not working DrawFunctions(functions_normal)\n"
                << "not working DrawFunctions(functions_ortog)\n";
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
#define EXECUTION_PART all

#if EXECUTION_PART == 10 || EXECUTION_PART == all
        DispAbsBeatwen__4X<double, double>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 800);
#endif
        ////Part No. 1: Quadrature formulas of Newton-Kot(e)ca and Gauss
        //1.1 calculations of a certain integral using compound quadrature formulas 
        configurate = example_2;
#if EXECUTION_PART == 11 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, IntegrateMethod::LEFT_RECTANGLE>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::MIDDLE_RECTANGLE>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::TRAPEZOID>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::SIMPSON>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
#endif

        //Methods for calculating a definite integral using compound quadrature formulas
        //based on 3 - point Newton - Kot formulas(e)ca and Gauss.
#if EXECUTION_PART == 12 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, IntegrateMethod::NEWTON_COTES_3_POINT> (configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::NEWTON_COTES_4_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::NEWTON_COTES_5_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
#endif
        //1.3. Graph of the dependence of the absolute error on the number of splits 
        //of the integration interval of each quadrature formula from paragraphs 1.1 - 2.
#if EXECUTION_PART == 13 || EXECUTION_PART == all
        DispAbsBeatwen<double, double, IntegrateMethod::GAUS_3_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 100);
        DispAbsBeatwen<double, double, IntegrateMethod::GAUS_4_POINT>(configurate.function, configurate.a, configurate.b, configurate.expected, { configurate.alpha,configurate.beta }, 20);
#endif
        ////Part No. 2: Methods for estimating the error of composite quadrature formulas
        //A definite integral with a given accuracy 𝜀 using the composite 3 - point quadrature formula of Newton - Cot(e)ca.
        // 
        //Estimation of the error by the Richardson method.
        //The rate of convergence according to Aitken's rule.
        // 
        //The length of the step ℎ of the partition.
        IntegrationResult<double> result;

        
#if EXECUTION_PART == 21 || EXECUTION_PART == all
        result = adaptive_integrate<IntegrateMethod::NEWTON_COTES_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, 1, { configurate.alpha, configurate.beta });
        std::cout <<"IntegrateMethod::NEWTON_COTES_3_POINT\n" << result;
#endif

        // N2.1 but using 3-point Gauss formulas:
#if EXECUTION_PART == 22 || EXECUTION_PART == all
        result = adaptive_integrate<IntegrateMethod::GAUS_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, 1, { configurate.alpha, configurate.beta });
        std::cout <<"IntegrateMethod::GAUS_3_POINT\n" << result;
#endif

        //Using the Aitken convergence rate estimate, select the h_opt step for either formulas 2.1 or 2.2.
        // 
        //Сompare it with the step calculated in 2.1 or 2.2.
#if EXECUTION_PART == 23 || EXECUTION_PART == all
        auto simple_result = optimized_adaptive_integration<IntegrateMethod::GAUS_3_POINT>(
                configurate.function, configurate.a, configurate.b, configurate.tolerance, { configurate.alpha, configurate.beta });
        std::cout << "IntegrateMethod::GAUS_3_POINT\n" << simple_result;
#endif
        ////
#undef EXECUTION_PART
    }



    namespace differential_equation_conf {

        template<typename Targ, typename Tresult>
        struct TestConfig {
            Targ c2;
            Targ x0;
            Tresult y0;
            Targ x_target;
            const uint64_t number_of_steps;

            std::function<Tresult(Targ, Tresult)> f;
            std::function<Tresult(Targ)> exact_y;

            double global_accuracy;
            double local_accuracy;
        };

        TestConfig<double, double> basic_scalar{
            0.5,

            0.0,
            1.0,

            M_PI,
            10,

            [](double x, double y) {return x * y; },
            [](double x) {return std::exp(x * x / 2); },
            1e-6,
            1e-6
        };
        
        TestConfig<double, matrix<double>> basic_vector{
    0.5,                    // c2
    0.0,                    // x0
    matrix<double>({{1.0}, {0.0}}),  // y0 - вектор-столбец 2x1 [y1, y2]^T
    M_PI,                   // x_target
    100,                    // number_of_steps

    // Функция f(x, y) возвращает вектор-столбец производных
    // Система: y1' = y2, y2' = -y1 (гармонический осциллятор)
    [](double x, matrix<double> y) -> matrix<double> {
                // y - вектор-столбец размером 2x1
                matrix<double> result(2, 1); // результат тоже вектор-столбец 2x1

                // y1' = y2
                result[0][0] = y[1][0];
                // y2' = -y1
                result[1][0] = -y[0][0];

                return result;
            },

            // Точное решение для системы y1' = y2, y2' = -y1
            // с начальными условиями y1(0)=1, y2(0)=0
            [](double x) -> matrix<double> {
                matrix<double> result(2, 1);
                result[0][0] = std::cos(x);  // y1(x) = cos(x)
                result[1][0] = -std::sin(x); // y2(x) = -sin(x)
                return result;
            },

            1e-6,   // global_accuracy
            1e-6    // local_accuracy
        };

        struct RightHandSideCounter {
            std::shared_ptr<size_t> call_count;

            RightHandSideCounter() : call_count(std::make_shared<size_t>(0)) {}

            matrix<double> operator()(double x, matrix<double> y) const {
                (*call_count)++;  

                matrix<double> result(2, 1);
                double A = 0.1;
                double B = 0.05;

                // dy1/dx = A * y2
                result[0][0] = A * y[1][0];
                // dy2/dx = -B * y1
                result[1][0] = -B * y[0][0];

                return result;
            }

            void reset() {
                *call_count = 0;
            }

            size_t count() const {
                return *call_count;
            }
        };

        RightHandSideCounter rhs_counter;
        TestConfig<double, matrix<double>> training_example{
            0.1,                    // c2 = 1/10
            0.0,                    // x0
            matrix<double>({{M_PI * 0.05}, {M_PI * 0.1}}),  // y0 = [B*π; A*π] = [π/20; π/10]
            M_PI,                   // x_target = π
            100,                    // number_of_steps

            // System: y1' = A*y2, y2' = -B*y1
            //  A = 1/10 = 0.1, B = 1/20 = 0.05
             std::ref(rhs_counter),

            [](double x) -> matrix<double> {
                double A = 0.1;
                double B = 0.05;
                double omega = std::sqrt(A * B); // ω = √(A*B)

                matrix<double> result(2, 1);

                // The exact solution of the systemy1' = A*y2, y2' = -B*y1
                // with initial conditions y1(0)=B*π, y2(0)=A*π
                result[0][0] = M_PI * B * std::cos(omega * x) +
                              M_PI * (A * A / omega) * std::sin(omega * x);

                result[1][0] = M_PI * A * std::cos(omega * x) -
                              M_PI * (B * omega / A) * std::sin(omega * x);

                return result;
            },

            1e-6,   // global_accuracy
            1e-6    // local_accuracy
        };
    }
    using namespace differential_equation_conf;

    template<typename DifferencialMethod method1, typename DifferencialMethod method2, typename Targ, typename Tresult>
    void CauchyStepSizeComparasion(TestConfig<Targ, Tresult> conf) {
        auto compute_error = [](const auto& var, const auto& working_example) {
            using ValueType = std::decay_t<decltype(var.second)>;

            if constexpr (std::is_arithmetic_v<ValueType>) {
                return var.second - working_example.exact_y(var.first);
            }
            else if constexpr (is_matrix<ValueType>::value) {
                return (var.second - working_example.exact_y(var.first)).norm();
            }
            };

        auto &working_example = conf;

        std::stringstream buffer;

        auto solution4 = SolveODE<method1, ErrorMethod::RungeGlobal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, working_example.global_accuracy, working_example.c2);

        buffer << "sol4_RK2 \n";

        for (auto& var : solution4)
        {
            //buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}" << "difference:" << compute_error(var, working_example) << "\n";
        }
        auto var = solution4.back();

        buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
            << "difference:" << compute_error(var, working_example) << "\n";

        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto h_eps = GetEpsStepSize(
            (solution4[2].first - solution4[1].first),
            1e-4,
            RungeError(solution4.back().second, (
                details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method1>(
                    working_example.f,
                    working_example.x0,
                    working_example.y0,
                    working_example.x_target,
                    solution4.size() * 2,
                    working_example.c2)
                ).back().second,
                order_of_accuracy_differencial_method.at(method1)),

            order_of_accuracy_differencial_method.at(method1));

        std::cout << "GetEpsStepSize " << h_eps << "\n";

        solution4 = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method1>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, std::ceil((working_example.x_target - working_example.x0) / h_eps), working_example.c2);

        std::vector<std::pair<double, double>> display_vec_1;
        buffer << "sol4_RK2 h_opt \n";
        for (auto& var : solution4)
        {
            //buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}" << "difference:" << compute_error(var, working_example) << "\n";

            display_vec_1.push_back({ var.first , compute_error(var, working_example) });
        }
        var = solution4.back();

        buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
            << "difference:" << compute_error(var, working_example) << "\n";

        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto getErrorIfLessOrEqual1 = [display_vec_1](double input_value) ->double {
            for (const auto& pair : display_vec_1) {
                if (input_value <= pair.first) {
                    return pair.second;
                }
            }
            return 0.0;
            };
        ////////////////
        solution4 = SolveODE<method2, ErrorMethod::RungeGlobal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, working_example.global_accuracy, working_example.c2);

        buffer << "\n\nsol4_RK3 :\n";
        for (auto& var : solution4)
        {
            //buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}" << "difference:" << compute_error(var, working_example) << "\n";
        }
        var = solution4.back();

        buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
            << "difference:" << compute_error(var, working_example) << "\n";


        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto h_eps_2 = GetEpsStepSize(
            (solution4[2].first - solution4[1].first),
            1e-5,
            RungeError(solution4.back().second, (
                details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method2>(
                    working_example.f,
                    working_example.x0,
                    working_example.y0,
                    working_example.x_target,
                    solution4.size() * 2,
                    working_example.c2)
                ).back().second,
                order_of_accuracy_differencial_method.at(method2)),

            order_of_accuracy_differencial_method.at(method2));

        std::cout << "GetEpsStepSize " << h_eps_2 << "\n";
        std::cout << "(working_example.x_target - working_example.x0)" << (working_example.x_target - working_example.x0) << "\n";
        std::cout << " std::ceil((working_example.x_target - working_example.x0) / h_eps_2)" <</* std::ceil*/((working_example.x_target - working_example.x0) / h_eps_2) << "\n";

        solution4 = SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method2>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, std::ceil((working_example.x_target - working_example.x0) / h_eps_2), working_example.c2);

        std::vector<std::pair<double, double>> display_vec_2;

        buffer << "sol4_RK3 h_opt :\n";

        std::cout << "\n\n\n\n\n\n\n\n mY solution4 size!!!!!!!!!!!!!!! " << solution4.size() << "\n";
        for (auto& var : solution4)
        {
            //buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}" << "difference:" << compute_error(var, working_example) << "\n";
            display_vec_2.push_back({ var.first ,compute_error(var, working_example) });
        }
        var = solution4.back();

        buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
            << "difference:" << compute_error(var, working_example) << "\n";

        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();


        auto getErrorIfLessOrEqual2 = [display_vec_2](double input_value) -> double {
            for (const auto& pair : display_vec_2) {
                if (0 <= input_value && input_value <= pair.first) {
                    return pair.second;
                }
            }
            return 0.0;
            };
        std::vector <std::function<double(double) >> functions;
        functions.push_back(getErrorIfLessOrEqual1);
        functions.push_back(getErrorIfLessOrEqual2);

        DrawFunctions(functions, 0.0, M_PI);
    }

    template<typename DifferencialMethod method1, typename DifferencialMethod method2, typename Targ, typename Tresult>
    void CauchyComparasionRatioTrueLocalErrorToEstimatedLocalError_vs_X(TestConfig<Targ, Tresult> conf)
    {
        std::stringstream buffer;
        auto& working_example = conf;
        auto compute_error = [](const auto& var, const auto& working_example) {
            using ValueType = std::decay_t<decltype(var.second)>;

            if constexpr (std::is_arithmetic_v<ValueType>) {
                return var.second - working_example.exact_y(var.first);
            }
            else if constexpr (is_matrix<ValueType>::value) {
                return (var.second - working_example.exact_y(var.first)).norm();
            }
            };


        auto solution6 = SolveODE<method1, ErrorMethod::RungeLocal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, 1e-9, working_example.c2);

        std::vector<std::pair<double, double>> vector_x_ratior_true_err_and_estimated_err_method1;
        for (int i = 1; i < solution6.size(); i++)
        {
            auto var = solution6[i];
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << solution6[i].first - solution6[i - 1].first << "\n";

            vector_x_ratior_true_err_and_estimated_err_method1.push_back({ solution6[i].first,compute_error(var, working_example) / (
                RungeError(details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method1>(working_example.f, solution6[i - 1].first,solution6[i - 1].second,solution6[i].first,2,working_example.c2).back().second,
                    details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method1>(working_example.f, solution6[i - 1].first,solution6[i - 1].second,solution6[i].first,4,working_example.c2).back().second,
                    order_of_accuracy_differencial_method.at(method1)
                    )) });
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto getOXhFunction1 = [vector_x_ratior_true_err_and_estimated_err_method1](double input_value) -> double {
            for (const auto& pair : vector_x_ratior_true_err_and_estimated_err_method1) {
                if (0 <= input_value && input_value <= pair.first) {
                    return pair.second;
                }
            }
            return 0.0;
            };
        ////
        solution6 = SolveODE<method2, ErrorMethod::RungeLocal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, 1e-9, working_example.c2);

        std::vector<std::pair<double, double>> vector_x_ratior_true_err_and_estimated_err_method2;
        for (int i = 1; i < solution6.size(); i++)
        {
            auto var = solution6[i];
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << solution6[i].first - solution6[i - 1].first << "\n";

            vector_x_ratior_true_err_and_estimated_err_method2.push_back({ solution6[i].first,compute_error(var, working_example) / (
                RungeError(details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method2>(working_example.f, solution6[i - 1].first,solution6[i - 1].second,solution6[i].first,2,working_example.c2).back().second,
                    details::SolveASystemOfOrdinaryDifferentialEquationsEqualSteps<method2>(working_example.f, solution6[i - 1].first,solution6[i - 1].second,solution6[i].first,4,working_example.c2).back().second,
                    order_of_accuracy_differencial_method.at(method2)
                    )) });
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto getOXhFunction2 = [vector_x_ratior_true_err_and_estimated_err_method2](double input_value) -> double {
            for (const auto& pair : vector_x_ratior_true_err_and_estimated_err_method2) {
                if (0 <= input_value && input_value <= pair.first) {
                    return pair.second;
                }
            }
            return 0.0;
            };


        ///
        std::vector<std::function<double(double)>> getOXhFunctionsVector;
        getOXhFunctionsVector.push_back({ getOXhFunction1 });
        getOXhFunctionsVector.push_back({ getOXhFunction2});
        DrawFunctions(getOXhFunctionsVector, 0.0, M_PI);


    }

    template<typename DifferencialMethod method1, typename DifferencialMethod method2, typename Targ, typename Tresult>
    void CauchyRatioTrueLocalAndEstimatedErrorsComparasion(TestConfig<Targ, Tresult> conf) {
        auto compute_error = [](const auto& var, const auto& working_example) {
            using ValueType = std::decay_t<decltype(var.second)>;

            if constexpr (std::is_arithmetic_v<ValueType>) {
                return var.second - working_example.exact_y(var.first);
            }
            else if constexpr (is_matrix<ValueType>::value) {
                return (var.second - working_example.exact_y(var.first)).norm();
            }
            };
    }


    void sem_5_part2(){
#define EXECUTION_PART 333
        /*
        The Cauchy problem :
            dy1(𝑥) / dx = 𝐴_𝑦2(𝑥),
            dy2(𝑥) / dx = -B_𝑦1(𝑥),
            y1(0) = Bπ,
            y2(0) = Aπ,
        */

        auto working_example = training_example;
        std::stringstream buffer;


        auto compute_error = [](const auto& var, const auto& working_example) {
            using ValueType = std::decay_t<decltype(var.second)>;

            if constexpr (std::is_arithmetic_v<ValueType>) {
                return var.second - working_example.exact_y(var.first);
            }
            else if constexpr (is_matrix<ValueType>::value) {
                return (var.second - working_example.exact_y(var.first)).norm();
            }
            };

        ////Part No. 1: Calculation schemes of the Runge-Kutta method with a constant step
        //A second-order calculation scheme based on second-order conditions for the 2-stage explicit Rungi-Kutta method for the c2 parameter
#if EXECUTION_PART == 11 || EXECUTION_PART == all
       
         

       


#endif  

        // Implementation of RK scheme construction with parameter c2 = ξ
        // 1.2 Implement constant-step RK method with total error estimation using Runge method
        //     (ε = 1e-4), initial step selection according to algorithm
#if EXECUTION_PART == 12 || EXECUTION_PART == all
        
        auto solution1 = SolveASystemOfOrdinaryDifferentialEquationsEqualStepsWithATotalError<DifferencialMethod::RungeKutta3ndOrder1,double,matrix<double>>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, working_example.global_accuracy, working_example.c2);
        std::cout << "sol1 \n";
        
        for (auto& var : solution1)
        {
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << compute_error(var, working_example) << "\n";
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();
        //std::cout << RungeError(solution.back().second, solution2.back().second,2);
#endif
        ////

        ////Part No. 2: Calculation schemes of the Runge-Kutta method with a constant step
        // 2.1 Implement automatic step-size control using 2-stage RK 2nd order method
        //     with local error estimation (ρ = 1e-5) and Runge method
#if EXECUTION_PART == 21 || EXECUTION_PART == all
        auto solution2 = SolveASystemOfOrdinaryDifferentialEquationsDifferentStepsWithALocalError< DifferencialMethod::RungeKutta2ndOrder>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, working_example.local_accuracy, working_example.c2);
        
        std::cout << "sol2 \n";
        for (auto& var : solution2)
        {
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << compute_error(var,working_example) << "\n";
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();
#endif
        ////

        ////Part No. 3: Calculation schemes of the Runge-Kutta method with a constant step
        // 3.1 Implement constant-step and automatic-step methods using classical RK schemes
        //     of 3rd or 4th order (opponent scheme)
#if EXECUTION_PART == 31 || EXECUTION_PART == all
        
        auto solution3 = SolveODE<DifferencialMethod::RK3, ErrorMethod::RungeGlobal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, working_example.local_accuracy, working_example.c2);

        std::cout << "sol3 \n";
        for (auto& var : solution3)
        {
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << compute_error(var, working_example) << "\n";
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

#endif

        // 3.2 Determine integration step h for constant-step methods (2-stage RK 2nd order
        //     and opponent scheme) that provides solution with accuracy ε = 1e-4
        //     Plot true total error vs x
#if EXECUTION_PART == 32 || EXECUTION_PART == all 
        CauchyStepSizeComparasion<DifferencialMethod::RK2, DifferencialMethod::RK3>(working_example);
#endif

        // 3.3 For automatic step-size control methods (2-stage RK 2nd order and opponent):
        // Analyze reliability and efficiency of implemented algorithms
        // 3.3.1 Plot step size vs x
#if EXECUTION_PART == 331 || EXECUTION_PART == all
        auto solution5 = SolveODE<DifferencialMethod::RK3, ErrorMethod::RungeLocal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, 1e-9, working_example.c2);
        
        std::cout << "sol331 \n" << solution5.size() << "\n";
        std::vector<std::pair<double, double>> oxh_display;
        for (int i = 0; i < solution5.size()-1;i++)
        {
            auto var = solution5[i];
            buffer << "{" << var.first << ";" << var.second << "} exact sol:{" << var.first << ";" << working_example.exact_y(var.first) << "}"
                << "difference:" << solution5[i + 1].first - solution5[i].first << "\n";
            oxh_display.push_back({ var.first+0.5,solution5[i + 1].first - solution5[i].first });
        }
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();

        auto getOXhFunction= [oxh_display](double input_value) -> double {
            for (const auto& pair : oxh_display) {
                if (0 <= input_value && input_value <= pair.first) {
                    return pair.second;
                }
            }
            return 0.0;
            };

        std::vector<std::function<double(double)>> getOXhFunctionsVector;
        getOXhFunctionsVector.push_back({getOXhFunction});
        DrawFunctions( getOXhFunctionsVector,0.0,M_PI );

#endif

        // 3.3.2 Plot ratio of true local error to estimated local error vs x
#if EXECUTION_PART == 332 || EXECUTION_PART == all

        CauchyComparasionRatioTrueLocalErrorToEstimatedLocalError_vs_X<DifferencialMethod::RK2, DifferencialMethod::RK3>(working_example);
        
#endif

        // 3.3.3 Plot number of right-hand side evaluations vs accuracy ε
#if EXECUTION_PART == 333 || EXECUTION_PART == all
        std::vector<double> display_call_count_number;
        for (int i=0; i < 15; i++) {

            auto solution = SolveODE<DifferencialMethod::RK3, ErrorMethod::RungeLocal>(working_example.f, working_example.x0, working_example.y0, working_example.x_target, 1.0/std::pow(2,3*i), working_example.c2);
            std::cout << "i:" << i << "\tcount:" << rhs_counter.count() << "\n";
            display_call_count_number.push_back(rhs_counter.count());
            rhs_counter.reset();

        }
        
        std::vector<std::function<double(double)>> getOXhFunctionsVector;
        getOXhFunctionsVector.push_back({ createLambda(display_call_count_number) });
        DrawFunctions(getOXhFunctionsVector);

#endif
        ////

#undef EXECUTION_PART

    }
}

/////////////////////////////////////////////////////////


int main() {

    

    //sem_4::sem_4();
    //sem_5::sem_5_part1();
    sem_5::sem_5_part2();



    system("pause");
    return 0;
}



