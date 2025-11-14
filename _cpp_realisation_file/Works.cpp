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
    counting_methods_2::polynomial_interpolation::nuton2::ShowInterpolationStatistic(func,n,m+50, #func, __VA_ARGS__)


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
        counting_methods_2::polynomial_interpolation::nuton2::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );


        counting_methods_2::aproximate::ShowAproximateStatistic([](double x) { return (x - std::sin(x) - 0.25); });

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
        counting_methods_2::polynomial_interpolation::nuton2::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<3, 2, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<3, 2, double>",
            [](double x) { return  identifier; },
            the_left_border, the_right_border,
            5
        );

        counting_methods_2::polynomial_interpolation::nuton2::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
            Spline_interpolator<2, 1, double>, // interpolator
            50, 50,                      // n, m_
            "Spline_interpolator<2, 1, double>",
            [](double x) { return identifier; },
            the_left_border, the_right_border,
            20
        );

        counting_methods_2::polynomial_interpolation::nuton2::ShowInterpolationStatistic<double, SplineInterpolatorFunc>(
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
        counting_methods_2::aproximate::ShowAproximateStatistic([](double x) { return (x - std::sin(x) - 0.25); });
#endif




    }
}

namespace sem_5 {
    void sem_5() {
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
//location_of_the_split_points in diapazon [0,1]
template<typename IntegratedFunction>
void GetaPartialIntegralSum_Singlethreaded(const IntegratedFunction integrated_function,
    const double left_boundary_of_the_integration_gap,
    const double right_boundary_of_the_integration_gap,
    uint64_t number_of_partitioning_intervals,
    matrix<double> location_of_the_split_points
) {
    if (location_of_the_split_points.is_vector() == false || (left_boundary_of_the_integration_gap > right_boundary_of_the_integration_gap)) { throw std::exception("interval_partitioning_scheme is not a vector "); }
    if (location_of_the_split_points.is_vertical_vector()) {
        location_of_the_split_points = location_of_the_split_points.transpose();
    }
    double ans_sum = 0;
    const double interval_size = (right_boundary_of_the_integration_gap - left_boundary_of_the_integration_gap) / number_of_partitioning_intervals;
    const double number_of_split_points = location_of_the_split_points.getrow();

    matrix<double> interval_partitioning_scheme = counting_methods_3::СoefficNewtonCotes(location_of_the_split_points, static_cast<double>(0), static_cast<double>(1));
    std::cout << "interval_partitioning_scheme" << interval_partitioning_scheme << '\n';



    ////nodes constructor
    std::vector<double> nodes_of_the_partition_on_the_interval(number_of_split_points, 0);
    auto F = [left_boundary_of_the_integration_gap, interval_size, number_of_split_points](double s) {return left_boundary_of_the_integration_gap + (interval_size * s); };

    /*nodes_of_the_partition_on_the_interval[0] = left_boundary_of_the_integration_gap;
    std::cout << "nodes_of_the_partition_on_the_interval[0]=" << nodes_of_the_partition_on_the_interval[0] << '\n';*/
    if (number_of_split_points > 1) {
        for (unsigned int i = 0; i < number_of_split_points; i++)
        {
            nodes_of_the_partition_on_the_interval[i] = F(location_of_the_split_points[0][i]);
            std::cout << "nodes_of_the_partition_on_the_interval[i]=" << nodes_of_the_partition_on_the_interval[i] << '\n';
        }
    }
    ////

    ////
    for (uint64_t i = 0; i < number_of_partitioning_intervals; i++) {
        double local_sum = 0;
        for (uint64_t j = 0; j < number_of_split_points; j++) {
            local_sum += integrated_function(nodes_of_the_partition_on_the_interval[j]) * interval_partitioning_scheme[0][j];

        }
        std::cout << "local_sum=" << local_sum << '\n';
        ans_sum += local_sum;
        for (auto i : nodes_of_the_partition_on_the_interval)
        {
            i += interval_size;
            std::cout << "nodes_of_the_partition_on_the_interval[j]=" << i << '\n';
        }
    }
    ////
    std::cout << "ans_sum=" << ans_sum << '\n';
}