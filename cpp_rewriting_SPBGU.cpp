#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "file_h/complex.h"
#include "file_h/fraction.h"
#include "file_h/polynomial.h"
#include "file_h/matrix.h"

#include <functional>
namespace calc_computing_f
{
void complex_calc() {

    static char L = 'E';
    while (1) {


        Complex a, b, c;
        std::cout << "Enter first Number \n";
        std::cin >> a;
        std::cout << "Enter second Number \n";
        std::cin >> b;

        std::cout << "Enter Action('-','+','/','*','E') :";

        std::cin >> L;
        if (L == 'E' || L == 'e')
        {
            break;
        }
        std::cout << std::endl;
        switch (L)
        {
        case '+':c = a + b; break;
        case '-':c = a - b; break;
        case '*':c = a * b; break;
        case '/':c = a / b; break;
        default:throw std::invalid_argument("unknown symbol");
        }
        std::cout << c << "\n";
    }


}
template<typename T>void matrix_calc() {

    static char L = 'E';
    while (1) {


        matrix<T> a, b, ans_out;
        T ans_out_T;
        std::cout << "Enter first Number \n";

        std::cin >> a;
        std::cout << "\n" << a << "\n";
        std::cout << "Enter second Number \n";
        std::cin >> b;
        std::cout << "\n" << b << "\n";
        std::cout << "Enter Action('-','+','/','*','r','E','d') :\n";
        std::cin >> L;
        if (L == 'E' || L == 'e')
        {
            break;
        }
        std::cout << std::endl;
        switch (L)
        {
        case '+':ans_out = a + b; break;
        case '-':ans_out = a - b; break;
        case '*':ans_out = a * b; break;
        case '/':ans_out = a / b; break;
        case 'd':ans_out_T = (a.determinant()); break;
        case 'r':ans_out = a.inverse_M(); break;
        default:throw std::invalid_argument("unknown symbol");
        }
        if (L == 'd') {
            std::cout << ans_out_T << "\n";
        }
        else
        {
            std::cout << ans_out << "\n";
        }


    }


}
void main_calc_menu() {
    std::cout << "Enter calc type \n";

    std::cout << "mfpint-matrix_calc<fraction<polynomial<int>>>()\n";
    std::cout << "mfpdouble-matrix_calc<fraction<polynomial<double>>>()\n";
    std::cout << "mpint-matrix_calc<polynomial<int>>()\n";
    std::cout << "mpdouble-matrix_calc<polynomial<double>>()\n";
    std::cout << "mdouble-matrix_calc<double>()\n";
    std::cout << "mint-matrix_calc<int>()\n";
    std::cout << "cmpl-complex_calc()\n";
    std::cout << "\ncalc type:";

    std::string type;
    std::cin >> type;

    if (type == "complex_calc()" || type == "cmpl")complex_calc();
    else if (type == "matrix_calc<fraction<polynomial<int>>>()" || type == "mfpint" || type == "matfrpolint") matrix_calc<fraction<polynomial<int>>>();
    else if (type == "matrix_calc<fraction<polynomial<double>>>()" || type == "mfpdouble" || type == "matfrpoldouble") matrix_calc<fraction<polynomial<double>>>();
    else if (type == "matrix_calc<polynomial<int>>()" || type == "mpint" || type == "matpolint") matrix_calc<polynomial<int>>();
    else if (type == "matrix_calc<polynomial<double>>()" || type == "mpdouble" || type == "matpoldouble") matrix_calc<polynomial<double>>();
    else if (type == "matrix_calc<double>()" || type == "mdouble" || type == "matdouble") matrix_calc<double>();
    else if (type == "matrix_calc<int>()" || type == "mint" || type == "matint") matrix_calc<int>();
    else throw std::invalid_argument("unknown calculator type or an error in the name");
}
//global computing function 

void calc_global() {
    try
    {
        std::string type;
    restart:
        main_calc_menu();
        std::cout << "Do you want to get out?(yes/no)\n:";
        std::cin >> type;
        if (type == "yes" || type == "y" || !(type == "no" || type == "n"))std::exit(0);
        if (type == "no" || type == "n")goto restart;

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
                std::cout<<next_guess[i]<<"\n";
            }

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

    int run_nnssi_with_setted_function() {
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

        void gaus_solver_linear_sistem() {
           

            try
            {
                matrix<double> a;
            std::cin >> a;
            std::cout << "\n" << a;

            matrix<double> b;
            std::cin >> b;
            std::cout << "\n" << b;

            std::cout << "\na*b:\n" << (a.inverse_M()) * b;

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
}


int main() {
    //calc_computing_f::matrix_calc<double>();


//C:\Users\User\source\repos\cpp_rewriting_SPBGU\input_matrix.txt
    /*counting_methods::polinomial::polynomial_test();*/
    //counting_methods::nonlinear_system_with_simple_iterations::run_nnssi_with_setted_function();

    //counting_methods::nonlinear_system_with_the_tangent_method::nonlinsystem_tangent_method();
    counting_methods::gaus_method::gaus_solver_linear_sistem();

    system("pause");
    return 0;
}
