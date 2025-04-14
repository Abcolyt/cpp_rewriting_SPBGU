#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <corecrt_math_defines.h>
#include <cmath>

#include "file_h/complex.h"
#include "file_h/fraction.h"
#include "file_h/polynomial.h"
#include "file_h/matrix.h"

#include "file_h/counting_methods_1.h"
#include "file_h/counting_methods_2.h"
#include "file_h/Array_xy_To_.h"
#include "file_h/calc_computing_f.h"





namespace counting_methods {
                                                          
    void executeWithFileInput(std::function<void()> func, const char* filename) {
        FILE* file;
        
        if (freopen_s(&file, filename, "r", stdin) != 0) {
            std::cerr << "Ошибка при открытии файла: " << filename << std::endl;
            return; // Завершение функции в случае ошибки
        }

        func();

        fclose(file);
    }




    namespace polinomial {


    template<typename P>P abs(P arg) {
        if (arg < 0) { arg = arg * (-1); }
        return arg;
    }

    template<typename P>std::pair<P, int> solve_tangents(polynomial<P> plnm,P x0, double e) {
        int iterations = 0;
        P x_old = 5, x_new = x0;
        while (abs(x_new - x_old) > e) {
            x_old = x_new;
            x_new = x_old +  (polynomialfunctions::f_polyn_x0_<P>(plnm, x_old) / (polynomialfunctions::f_polyn_x0_<P>(polynomialfunctions::derivate(plnm), x_old)))*(-1);
            iterations++;
        }

        return { x_new, iterations };
    }

    template<typename P>std::vector<std::pair<P, int>> plnm_roots(polynomial<P> plnm,P x0) {
        std::vector<std::pair<P,int>> ans_roots;
        polynomial<P> b;
        size_t deg=plnm.get_deg()-1;
        for (size_t i = 0; i < deg; i++)
        {
            auto ans = solve_tangents<P>(plnm, x0, LDBL_EPSILON);
            //std::cout <<"root["<<i<<"]=" << ans.first << '\n';
            b.newsize(2);
            b[1] = 1;
            b[0] = (ans.first) * (-1);

            plnm = (plnm / b);
            //std::cout << "plnm:" << plnm << "\nnew plnm" << b << "\nnew iter" << "'"<< i<< "' \n";
            ans_roots.push_back(ans);
        }
    


        return ans_roots;
    }


    void polynomial_test() {
        polynomial<Complex> plnm, b;
        std::cin >> plnm;
        auto ans = plnm_roots<Complex>(plnm,Complex(LDBL_EPSILON,LDBL_EPSILON));
        for (size_t i = 0; i < ans.size(); i++)
        {
            std::cout << "root:(" << ans[i].first << ")\niteration for this root:" << ans[i].second << " \n";
            //std::cout << "polynomialfunctions::f_polyn_x0_<Complex>(plnm, ans["<<i<<"].first):" << polynomialfunctions::f_polyn_x0_<Complex>(plnm, ans[i].first)<<'\n';
        }


    }


    }

    namespace nonlinear_system_with_simple_iterations {
    using Function = std::function<double(const std::vector<double>&)>;

    // Функция для решения нелинейной системы уравнений методом простых итераций
    std::vector<double> nonlinear_system_with_simple_iterations(
        const std::vector<Function>&functions,
        const std::vector<double>&initial_guess,
        double tolerance = 1e-20,
        int max_iterations = LONG_MAX)
    {
        std::vector<double> current_guess = initial_guess;
        int num_functions = functions.size();

        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            std::vector<double> next_guess(num_functions);

            // Вычисляем значения для следующей итерации
            for (int i = 0; i < num_functions; ++i) {
                next_guess[i] = functions[i](current_guess);
                std::cout<<next_guess[i]<<"    ";
            }
            std::cout << "\n";
            // Проверяем на сходимость
            double max_diff = 0.0;
            for (int i = 0; i < num_functions; ++i) {
                max_diff = std::max(max_diff, std::abs(next_guess[i] - current_guess[i]));
            }

            if (max_diff < tolerance) {
                return next_guess; // Возвращаем найденное решение
            }

            current_guess = next_guess; // Переходим к следующей итерации
        }

        throw std::runtime_error("Maximum iterations reached without convergence");
    }

    int run_nnssi_with_setted_nonlinear_function() {
        //We define the functions of the system of equations 13
        std::vector<Function> functions1 = {
            [](const std::vector<double>& x) { return 1 - ((1+LDBL_EPSILON) / (2)) * std::sin(x[1] + 1); },//x=1-((1)/(2))sin(y+1)
            [](const std::vector<double>& x) { return 0.7 - (1-LDBL_EPSILON)*cos(x[0] - 1); }                   //y=0.7-cos(x-1)
        };

        std::vector<double> initial_guess = { 0, 0};

        try {
            
            std::vector<double> solution = nonlinear_system_with_simple_iterations(functions1, initial_guess);

            //  результат
            std::cout << "Solution: ";
            for (double value : solution) {
                std::cout << value << " ";
            }
            std::cout << std::endl;

        }
        catch (const std::exception& e) {
            std::cerr << "!Error: " << e.what() << std::endl;
        }

        return 0;   
    }
    int run_nnssi_with_setted_linear_function() {
        //We define the functions of the system of equations 13
        std::vector<Function> functions1 = {
            [](const std::vector<double>& x) { return (24.4781 -(0.0496 * x[1] + 0.0444 * x[2] + 0.0393 * x[3])) / 16;},
            [](const std::vector<double>& x) { return (26.0849 -(.0688 * x[0] + .0585 * x[2] + .0534 * x[3])) / 15.1; } ,
            [](const std::vector<double>& x) { return (27.3281 -(.0829 * x[0] + .0777 * x[1] + .0674 * x[3])) / 14.2000; },
            [](const std::vector<double>& x) { return (28.2078 -(.0970 * x[0] + .0918 * x[1] + .0867 * x[2])) / 13.3000;  }
        };

        std::vector<double> initial_guess = { 0.1, 0.1,0.1,0.1 };

        try {

            std::vector<double> solution = nonlinear_system_with_simple_iterations(functions1, initial_guess);

            //  результат
            std::cout << "Solution: ";
            for (double value : solution) {
                std::cout << value << " ";
            }
            std::cout << std::endl;

        }
        catch (const std::exception& e) {
            std::cerr << "!Error: " << e.what() << std::endl;
        }

        return 0;
    }


    }





    namespace nonlinear_system_with_the_tangent_method {

        // Определим типы для функций и якобиана
        using Function = std::function<double(const std::vector<double>&)>;
        using Jacobian = std::function<std::vector<std::vector<double>>(const std::vector<double>&)>;

        // Функция для решения нелинейной системы уравнений методом Ньютона
        std::vector<double> nonlinear_system_with_newton(
            const std::vector<Function>& functions,
            const Jacobian& jacobian,
            const std::vector<double>& initial_guess,
            double tolerance = 1e-9,
            int max_iterations = 100000)
        {
            std::vector<double> current_guess = initial_guess;
            int num_functions = functions.size();

            for (int iteration = 0; iteration < max_iterations; ++iteration) {
                // Вычисляем значения функций
                std::vector<double> function_values(num_functions);
                for (int i = 0; i < num_functions; ++i) {
                    function_values[i] = functions[i](current_guess);
                }

                // Вычисляем якобиан
                std::vector<std::vector<double>> J = jacobian(current_guess);

                // Решаем систему J * delta = -F для delta
                std::vector<double> delta(num_functions);

                // Используем метод Гаусса или любой другой метод для решения системы
                // Для простоты, будем использовать метод Гаусса с прямой подстановкой
                // Здесь мы просто создаем матрицу и вектор для решения
                std::vector<std::vector<double>> augmented_matrix(num_functions, std::vector<double>(num_functions + 1));

                for (int i = 0; i < num_functions; ++i) {
                    for (int j = 0; j < num_functions; ++j) {
                        augmented_matrix[i][j] = J[i][j];
                    }
                    augmented_matrix[i][num_functions] = -function_values[i];
                }

                // Прямой ход Гаусса
                for (int i = 0; i < num_functions; ++i) {
                    // Нормализация строки
                    double pivot = augmented_matrix[i][i];
                    for (int j = i; j <= num_functions; ++j) {
                        augmented_matrix[i][j] /= pivot;
                    }

                    // Обнуление ниже
                    for (int k = i + 1; k < num_functions; ++k) {
                        double factor = augmented_matrix[k][i];
                        for (int j = i; j <= num_functions; ++j) {
                            augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
                        }
                    }
                }

                // Обратный ход Гаусса
                for (int i = num_functions - 1; i >= 0; --i) {
                    delta[i] = augmented_matrix[i][num_functions];
                    for (int j = i + 1; j < num_functions; ++j) {
                        delta[i] -= augmented_matrix[i][j] * delta[j];
                    }
                }

                // Обновляем текущее приближение
                for (int i = 0; i < num_functions; ++i) {
                    current_guess[i] += delta[i];
                }

                // Проверяем на сходимость
                double max_diff = 0.0;
                for (const auto& d : delta) {
                    max_diff = std::max(max_diff, std::abs(d));
                }

                if (max_diff < tolerance) {
                    return current_guess; // Возвращаем найденное решение
                }
            }

            throw std::runtime_error("Maximum iterations reached without convergence");
        }

        // Пример использования
        int nonlinsystem_tangent_method() {
            // Определяем функции системы уравнений
            std::vector<Function> functions = {
                [](const std::vector<double>& x) { return std::tan(x[0] * x[1] + 0.4) - x[0] * x[0]; }, // Пример: f1(x1, x2) = x1^2 + x2 - 2
                [](const std::vector<double>& x) { return 0.8 * x[0] * x[0] + 2 * x[1] * x[1] - 1; } // Пример: f2(x1, x2) = x1 - x2^2
            };

            // Определяем якобиан системы уравнений
            Jacobian jacobian = [](const std::vector<double>& x) {
                return std::vector<std::vector<double>>{
                    { x[1]/(std::cos((5*x[1]*x[0]+2)/5)* std::cos((5 * x[1] * x[0] + 2) / 5)) - 2*x[0] , x[0]/ (std::cos((5 * x[1] * x[0] + 2) / 5) * std::cos((5 * x[1] * x[0] + 2) / 5)) },     // df1/dx1, df1/dx2
                    {8*x[0],4*x[1]}     // df2/dx1, df2/dx2
                };
                };

            // Начальное приближение
            std::vector<double> initial_guess = {0.1, 0.1};

            try {
                // Решаем систему уравнений
                std::vector<double> solution = nonlinear_system_with_newton(functions, jacobian, initial_guess);

                // Выводим результат
                std::cout << "Solution: \n";
                for (double value : solution) {
                    std::cout << value << " presizion:   "<< std::sqrt(functions[0](solution) * functions[0](solution) + functions[1](solution)* functions[1](solution))<<"\n";
                }
                std::cout << std::endl;

            }
            catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << std::endl;
            }

            return 0;
        }
    }

    namespace gaus_method{

        matrix<double> mult(const matrix<double>& left, const matrix<double>& b)  {
            //if (left.getcol() != b.getrow() || b.getcol() != 1) {
            //    throw std::invalid_argument("Invalid dimensions for multiplication.");
            //}
            
            matrix<double> result(left.getrow(), 1);
            for (size_t i = 0; i < left.getrow(); ++i) {
                for (size_t j = 0; j < left.getrow(); ++j) {
                    result[0][i] += left[i][j] * b[0][j];
                    std::cout << "result[0][<<" << i << "]=" << result[0][i] << "\n";
                }
            }
            return result;
        }

        void gaus_solver_linear_sistem() {
           

            try
            {
                matrix<double> a;
            std::cin >> a;
            std::cout << "\n" << a;

            matrix<double> b;
            std::cin >> b;
            std::cout << "\n" << b;
            std::cout << "\na^-1:\n" << (a.inverse_M()) ;
            std::cout << "\na*b:\n" << mult((a.inverse_M()) , b);
            
            std::cout << "\ncheck:\n" << a * b;
            //std::cout << "\n(a.inverse_M()):\n" << (a.inverse_M());

            //std::cout << "\n a*(a.inverse_M()):\n" <<(a.inverse_M())* a;
            }



            catch (std::exception ex)
            {
                std::cout << "exeption!!What:" << ex.what();
            }
            catch (...)
            {
                std::cout << "unknown error";
            }
    }
    }
    
    namespace holechi {
        int example() {
        matrix<double> A(4);

        A[0][0] = .1954; A[0][1] = .7700; A[0][2] = 1.3446; A[0][3]= 1.9192;
        A[1][0] = .7700; A[1][1] = 15.1728; A[1][2] = 21.9666; A[1][3]= 28.7604;
        A[2][0] = 1.3446; A[2][1] = 21.9666; A[2][2] = 74.7291; A[2][3]= 93.3867;
        A[3][0] = 1.9192; A[3][1] = 28.7604; A[3][2] = 93.3867; A[3][3]= 208.6609;
        
        

            
        std::cout << "inpyt matrix:\n" << A;

        try {
            matrix<double> A_inv = (A.cholesky().inverse_M()).transpose()* A.cholesky().inverse_M();
            std::cout << " A^{-1}:\n" << A_inv;

            std::cout << "check A^{-1}:\n" << A_inv*A;
        }
        catch (const std::exception& e) {
            std::cerr << "error: " << e.what() << std::endl;
        }

        return 0;
    }
    }
}

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
    counting_methods_2::Polynomial_interpolation::nuton2::show_interpolation_statistic(func,n,m, #func, __VA_ARGS__)
template <typename Func, typename... Args>
constexpr void show_wrapper(Func func,double n, double m, Args... args) {
    SHOW_INTERPOL_STAT(func,n,m, std::forward<Args>(args)...);
}

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
#define AU_LOG 6
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
    /*show_interpolation_statistic(nuton_interpolation,10, 10,
        [](double x) { return std::cos(x) / std::sin(x) + x * x; },
        -2 * M_PI, 2 * M_PI
    );*/

SHOW_INTERPOL_STAT(
    nuton_interpolation<double>, // interpolator
    10, 20,                      // n, m_
    [](double x) { return std::cos(x) / std::sin(x) + x * x; },
    -2 * M_PI, 2 * M_PI
);
SHOW_INTERPOL_STAT(
    Lagrang_interpolation<double>, 
    10, 10,                      
    [](double x) { return std::cos(x) / std::sin(x) + x * x; },
    -2 * M_PI, 2 * M_PI
);
#endif
    
    system("pause");
    return 0;
}
