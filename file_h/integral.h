#pragma once
#include "../file_h/matrix.h"
#include <vector>
namespace counting_methods_3 {
#if 0
    template<typename T>
	matrix<T> СoefficNewtonCotes(std::vector<T> nodes_in_the_integration_gap, T beginning_of_the_integration_interval, T end_of_the_integration_interval) {
        
        matrix<T>A(nodes_in_the_integration_gap.size()), b(nodes_in_the_integration_gap.size(), 1);

        for (int i = 0; i < nodes_in_the_integration_gap.size(); i++) {
            A[0][i] = 1;//nodes_in_the_integration_gap[i];
            b[i][0] = 1 / (static_cast<T>(1 + i));
            for (int j = 1; j < nodes_in_the_integration_gap.size(); j++) {
                A[j][i] = nodes_in_the_integration_gap[i] * A[j - 1][i];
            }
        }
		
		matrix<T> coeff = (end_of_the_integration_interval - beginning_of_the_integration_interval) * matrixfunction::solve_system(A,b);

		return  coeff;
 }
#endif
    matrix<double> getrange(int n) {
        matrix<double> M;
        M.setcol(1);
        M.setrow(n);
        for (int i = 0; i < n; i++) {
            M[0][i] = i * (1.0 / (n - 1));
        }

        return M;
    }

    //get c[j] : Integral(a,b) of function =~= Sum(j=1,n){c[j]*F(x_j)} 
    template<typename T>
    matrix<T> СoefficIntegralNewtonCotes(int n) {
        matrix<T> nodes_in_the_integration_gap = getrange(n);

        if (nodes_in_the_integration_gap.is_vector()) {

            //std::cout << "nodes_in_the_integration_gap:\n" << nodes_in_the_integration_gap << "\n";
            if (nodes_in_the_integration_gap.getcol() == 1) { nodes_in_the_integration_gap = nodes_in_the_integration_gap.transpose(); }
            T beginning_of_the_integration_interval = 0;
            T end_of_the_integration_interval       = 1;


            //std::cout << "1:" << nodes_in_the_integration_gap <<'\n';
            if (beginning_of_the_integration_interval != 0 || end_of_the_integration_interval != 1) {
                //std::cout << "2:" << (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow()))) << "\n";
                nodes_in_the_integration_gap = (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow())));
                nodes_in_the_integration_gap = nodes_in_the_integration_gap /(end_of_the_integration_interval - beginning_of_the_integration_interval);
            }
            uint64_t size= std::max(nodes_in_the_integration_gap.getrow(), nodes_in_the_integration_gap.getcol());
            matrix<T>A(size), b(size, 1);
            //std::cout << "3:" << nodes_in_the_integration_gap << '\n';
            for (int i = 0; i < size; i++) {
                A[0][i] = 1;//nodes_in_the_integration_gap[i];
                b[i][0] = 1 / (static_cast<T>(1 + i));
                for (int j = 1; j < size; j++) {
                    A[j][i] = nodes_in_the_integration_gap[i][0] * A[j - 1][i];
                }
            }

            std::cout << "A:" << A << "\n";
            std::cout << "b:" << b << "\n";
            matrix<T> coeff = matrixfunction::solve_system(A, b);
            //std::cout << "с_i++:" << A.inverse_M()*b << "\n";
            
            std::cout << "c_i:" << coeff << '\n';
            return coeff;
        }
    }

    //get c[j] : Integral(a,b) of function =~= Sum(j=1,n){(b-a)*c[j]*F(x_j)} 
    template<typename T>
    matrix<T> СoefficNewtonCotes(matrix<T> nodes_in_the_integration_gap, T beginning_of_the_integration_interval, T end_of_the_integration_interval) {
        if (nodes_in_the_integration_gap.is_vector()) {
            //std::cout << "nodes_in_the_integration_gap:\n" << nodes_in_the_integration_gap << "\n";
            if (nodes_in_the_integration_gap.getcol() == 1) { nodes_in_the_integration_gap = nodes_in_the_integration_gap.transpose(); }
            std::cout <<"1:" << nodes_in_the_integration_gap << '\n';
            if (beginning_of_the_integration_interval != 0 || end_of_the_integration_interval != 1) {
                std::cout << "2:" << (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow()))) << "\n";
                nodes_in_the_integration_gap = (nodes_in_the_integration_gap - (matrix<T>::ones(nodes_in_the_integration_gap.getcol(), nodes_in_the_integration_gap.getrow())));
                nodes_in_the_integration_gap = nodes_in_the_integration_gap / (end_of_the_integration_interval - beginning_of_the_integration_interval);
            }
            uint64_t size = std::max(nodes_in_the_integration_gap.getrow(), nodes_in_the_integration_gap.getcol());
            matrix<T>A(size), b(size, 1);

            std::cout << "3:" << nodes_in_the_integration_gap << '\n';

            for (int i = 0; i < size; i++) {
                A[0][i] = 1;//nodes_in_the_integration_gap[i];
                b[i][0] = 1 / (static_cast<T>(1 + i));
                for (int j = 1; j < size; j++) {
                    A[j][i] = nodes_in_the_integration_gap[i][0] * A[j - 1][i];
                }
            }
            std::cout << "A:" << A<<"\n";
            std::cout << "b:" << b << "\n";

            matrix<T> coeff =  matrixfunction::solve_system(A, b);
            std::cout << "c_i:" << coeff << '\n';
            return coeff;
        }
    }


    
    ////
    // Calculating the logarithm of the beta function: ln(B(p, q))
    long double log_beta_function(long double p, long double q) {
        return std::lgamma(p) + std::lgamma(q) - std::lgamma(p + q);
    }
    //Calculating the beta function
    double beta_function(double p, double q) {
        return std::tgamma(p) * std::tgamma(q) / std::tgamma(p + q);
    }
    ////
    


    // Calculation of the integral I = ∫[a,b] x^j / ((x-a)^α * (b-x)^β) dx
    double compute_integral(double a, double b, double alpha, double beta, int j) {
        // Проверка корректности параметров
        if (a >= b) {
            throw std::invalid_argument("a must be less than b");
        }
        if (j < 0) {
            throw std::invalid_argument("j must be non-negative");
        }

        // Особый случай: a = 0
        if (a == 0.0) {
            double C1 = std::pow(b - a, 1 - alpha - beta + j);
            return C1 * beta_function(j - alpha + 1, 1 - beta);
        }

        // Общий случай: a ≠ 0
        double C1 = std::pow(b - a, 1 - alpha - beta);
        double C2 = std::pow(a, j);
        double C3 = (b - a) / a;

        double sum = 0.0;
        double binom_coeff = 1.0;  // Начинаем с C(j,0) = 1
        double C3_power = 1.0;     // Начинаем с C3^0 = 1

        for (int k = 0; k <= j; ++k) {
            // Вычисляем бета-функцию для текущего k
            double beta_val = beta_function(k - alpha + 1, 1 - beta);

            // Добавляем текущий член к сумме
            sum += binom_coeff * C3_power * beta_val;

            // Обновляем биномиальный коэффициент для следующей итерации
            if (k < j) {
                binom_coeff *= static_cast<double>(j - k) / (k + 1);
            }

            // Обновляем степень C3 для следующей итерации
            C3_power *= C3;
        }

        return C1 * C2 * sum;
    }

    // Calculating the integral using logarithmic calculations : 
    long double compute_integral_log(long double a, long double b,
        long double alpha, long double beta, int j) {
        // Проверка корректности параметров
        if (a >= b) {
            throw std::invalid_argument("a must be less than b");
        }
        if (j < 0) {
            throw std::invalid_argument("j must be non-negative");
        }

        // Особый случай: a = 0
        if (a == 0.0L) {
            long double log_C1 = (1 - alpha - beta + j) * std::log(b - a);
            long double log_beta = log_beta_function(j - alpha + 1, 1 - beta);
            return std::exp(log_C1 + log_beta);
        }

        // Общий случай: a ≠ 0
        long double log_C1 = (1 - alpha - beta) * std::log(b - a);
        long double log_C2 = j * std::log(a);
        long double log_C3 = std::log(b - a) - std::log(a);
        //


        // Вектор для хранения логарифмов членов суммы
        std::vector<long double> log_terms(j + 1);

        // Вычисляем логарифмы всех членов суммы
        long double log_binom = 0.0L; // ln(C(j,0)) = ln(1) = 0

        for (int k = 0; k <= j; ++k) {
            // Логарифм биномиального коэффициента + k*ln(C3)
            long double log_binom_C3 = log_binom + k * log_C3;

            // Логарифм бета-функции
            long double log_beta_val = log_beta_function(k - alpha + 1, 1 - beta);

            // Суммируем: ln(term_k) = ln(binom) + k*ln(C3) + ln(B(...))
            log_terms[k] = log_binom_C3 + log_beta_val;

            // Обновляем логарифм биномиального коэффициента для следующей итерации
            if (k < j) {
                // ln(C(j,k+1)) = ln(C(j,k)) + ln(j-k) - ln(k+1)
                log_binom += std::log(j - k) - std::log(k + 1);
            }
        }

        // Находим максимальный логарифм для численной стабильности
        long double max_log = log_terms[0];
        for (int k = 1; k <= j; ++k) {
            if (log_terms[k] > max_log) {
                max_log = log_terms[k];
            }
        }

        // Суммируем exp(log_terms[k] - max_log) и затем умножаем на exp(max_log)
        long double sum_exp = 0.0L;
        for (int k = 0; k <= j; ++k) {
            sum_exp += std::exp(log_terms[k] - max_log);
        }

        // Итоговый результат: C1 * C2 * sum
        long double log_result = log_C1 + log_C2 + max_log + std::log(sum_exp);
        return std::exp(log_result);
    }
    ////
    
    //Calculation of the integral I = ∫[a, b] x ^ j / ((x - a) ^ α * (b - x) ^ β) dx
    long double compute_integral_optimized(long double a, long double b,
        long double alpha, long double beta, int j) {



        if (a >= b) throw std::invalid_argument("a must be less than b");
        if (j < 0) throw std::invalid_argument("j must be non-negative");

        // Особый случай: a = 0
        if (a == 0.0L) {
            long double C1 = std::pow(b - a, 1 - alpha - beta + j);
            return C1 * std::exp(log_beta_function(j - alpha + 1, 1 - beta));
        }

        if (j <= 30) {
            // Прямой метод для малых j
            long double C1 = std::pow(b - a, 1 - alpha - beta);
            long double C2 = std::pow(a, j);
            long double C3 = (b - a) / a;

            long double sum = 0.0L;
            long double binom_coeff = 1.0L;
            long double C3_power = 1.0L;

            for (int k = 0; k <= j; ++k) {
                long double beta_val = std::exp(log_beta_function(k - alpha + 1, 1 - beta));
                sum += binom_coeff * C3_power * beta_val;

                if (k < j) {
                    binom_coeff *= static_cast<long double>(j - k) / (k + 1);
                }
                C3_power *= C3;
            }

            return C1 * C2 * sum;
        }
        else {
            // Общий случай: a ≠ 0
            long double log_C1 = (1 - alpha - beta) * std::log(b - a);
            long double log_C2 = j * std::log(a);
            long double log_C3 = std::log(b - a) - std::log(a);

            // Вектор для хранения логарифмов членов суммы
            std::vector<long double> log_terms(j + 1);

            // Вычисляем логарифмы всех членов суммы
            long double log_binom = 0.0L; // ln(C(j,0)) = ln(1) = 0

            for (int k = 0; k <= j; ++k) {
                // Логарифм биномиального коэффициента + k*ln(C3)
                long double log_binom_C3 = log_binom + k * log_C3;

                // log beta function
                long double log_beta_val = log_beta_function(k - alpha + 1, 1 - beta);

                // summ: ln(term_k) = ln(binom) + k*ln(C3) + ln(B(...))
                log_terms[k] = log_binom_C3 + log_beta_val;

                // Обновляем логарифм биномиального коэффициента для следующей итерации
                if (k < j) {
                    // ln(C(j,k+1)) = ln(C(j,k)) + ln(j-k) - ln(k+1)
                    log_binom += std::log(j - k) - std::log(k + 1);
                }
            }

            // Находим максимальный логарифм для численной стабильности
            long double max_log = log_terms[0];
            for (int k = 1; k <= j; ++k) {
                if (log_terms[k] > max_log) {
                    max_log = log_terms[k];
                }
            }

            // Суммируем exp(log_terms[k] - max_log) и затем умножаем на exp(max_log)
            long double sum_exp = 0.0L;
            for (int k = 0; k <= j; ++k) {
                sum_exp += std::exp(log_terms[k] - max_log);
            }

            // Итоговый результат: C1 * C2 * sum
            long double log_result = log_C1 + log_C2 + max_log + std::log(sum_exp);
            return std::exp(log_result);
        }
    }

    void test_integral_p_x_() {
        try {
            long double a = 1.0L, b = 3.0L;
            long double alpha = 0.5L, beta = 0.5L;

            std::cout << "Comparison of calculation methods:" << std::endl;
            std::cout << "a=" << a << ", b=" << b << ", α=" << alpha << ", β=" << beta << std::endl;

            for (int j = 0; j <= 10; j++) {
                long double result_log = counting_methods_3::compute_integral_log(a, b, alpha, beta, j);
                long double result_direct = counting_methods_3::compute_integral(a, b, alpha, beta, j);
                long double result_optim = counting_methods_3::compute_integral_optimized(a, b, alpha, beta, j);
                std::cout << "j=" << j << ": log_method=" << result_log
                    << ", direct_method=" << result_direct
                    << ", diff(log - opt)=" << std::abs(result_log - result_optim)
                    << ", diff(dir-opt)=" << std::abs(result_direct - result_optim) << std::endl;
            }

            // Тест с большими значениями j
            std::cout << "\nThe Big j test:" << std::endl;
            int large_j = 70;
            long double result_large = counting_methods_3::compute_integral_log(a, b, alpha, beta, large_j);
            long double opt_large = counting_methods_3::compute_integral_optimized(a, b, alpha, beta, large_j);
            std::cout << "j=" << large_j << ": result=" << result_large
                << ", diff+large(dir-opt)=" << std::abs(result_large - opt_large) << std::endl;

        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }

    }

    std::unordered_map<int, long double> fill_integral_map(
        long double a, long double b,
        long double alpha, long double beta,
        int n) {

        std::unordered_map<int, long double> integral_map;

        // for j = 0,..., 2*n-1
        for (int j = 0; j < 2 * n; ++j) {
            try {
                long double value = compute_integral_optimized(a, b, alpha, beta, j);
                integral_map[j] = value;

                std::cout << "j = " << j << ", I = " << value << std::endl;

            }
            catch (const std::exception& e) {
                std::cerr << "Error computing integral for j=" << j
                    << ": " << e.what() << std::endl;
                integral_map[j] = 0.0L; 
            }
        }

        return integral_map;
    }
    void fill_umap_test() {
        try {
            // Фиксированные параметры
            long double a = 1.0L;
            long double b = 3.0L;
            long double alpha = 0.5L;
            long double beta = 0.5L;
            int n = 5; // Будем вычислять для j = 0, 1, ..., 9 (2*5-1 = 9)

            std::cout << "Parameters: a=" << a << ", b=" << b
                << ", alpha=" << alpha << ", beta=" << beta
                << ", n=" << n << std::endl;
            std::cout << "Computing integrals for j = 0 to " << 2 * n - 1 << std::endl;

            // Заполняем unordered_map
            auto integral_map = counting_methods_3::fill_integral_map(a, b, alpha, beta, n);

            // Выводим результаты
            std::cout << "\nResults stored in unordered_map:" << std::endl;
            for (const auto& [j, value] : integral_map) {
                std::cout << "integral_map[" << j << "] = " << value << std::endl;
            }

            // Демонстрация доступа к элементам
            std::cout << "\nAccessing specific elements:" << std::endl;
            for (int j = 0; j < 2 * n; ++j) {
                if (integral_map.find(j) != integral_map.end()) {
                    std::cout << "integral_map[" << j << "] = " << integral_map[j] << std::endl;
                }
            }

        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            
        }
    }
    enum class IntegrateMethod {    
        LEFT_RECTANGLE,
        MIDDLE_RECTANGLE,
        TRAPEZOID,
        SIMPSON,
        NEWTON_COTES_3_POINT=3,
        NEWTON_COTES_4_POINT=4,
        NEWTON_COTES_5_POINT,
        NEWTON_COTES_6_POINT,
        NEWTON_COTES_7_POINT,
        NEWTON_COTES_8_POINT,
        NEWTON_COTES_9_POINT,
        GAUS_3_POINT,
        GAUS_4_POINT
    };
    // Исправленная сигнатура функции интегрирования
    template<typename TResult, typename TArg, IntegrateMethod Method>
    TResult integrate(std::function<TResult(TArg)> f, TArg a, TArg b, int points) {
        TResult sum = 0;
        if constexpr (Method == IntegrateMethod::LEFT_RECTANGLE) {
            TResult h = (b - a) / static_cast<TArg>(points);
            sum = 0;
            for (int i = 0; i < points; ++i) {
                sum += f(a + i * h);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::MIDDLE_RECTANGLE) {
            TArg h = (b - a) / static_cast<TArg>(points);
            sum = 0;
            for (int i = 1; i < points; ++i) {
                sum += f(a + i * h - h / 2);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::TRAPEZOID) {
            TArg h = (b - a) / static_cast<TArg>(points);
            sum = (f(a) + f(b)) / 2;
            for (int i = 1; i < points; ++i) {
                sum += f(a + i * h);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::SIMPSON) {
            TArg h = (b - a) / static_cast<TArg>(points);
            sum = f(a) + f(b);

            for (int i = 1; i < points; ++i) {
                TResult coefficient = (i % 2 == 1) ? 4 : 2;
                sum += coefficient * f(a + i * h);
            }

            return sum * h / 3;
        }
        else if constexpr ((static_cast<int>(IntegrateMethod::NEWTON_COTES_4_POINT) <= static_cast<int>(Method)) && (static_cast<int>(Method) <= static_cast<int>(IntegrateMethod::NEWTON_COTES_9_POINT))) {
            std::cout << "good" << '\n';
            TResult h = (b - a) / static_cast<TArg>(points);
            sum = 0;
            int n = static_cast<int>(Method);
            TResult micro_h = h / static_cast<TArg>(n-1);
            matrix<TResult> M = СoefficIntegralNewtonCotes<TResult>(n);

            std::cout<<"M:" << M;
            for (int j = 0; j < points; j++)
            {
                for (int i = 0; i < n; ++i) {
                   
                    sum += M[0][i] * f(a + h * j + micro_h * i);
                }
            }

            return sum * h;
        }
        //else if constexpr (Method == IntegrateMethod::GAUS_3_POINT) {

        //}
        //else if constexpr (Method == IntegrateMethod::GAUS_4_POINT) 

        //}

        // ... другие методы
        else {
            static_assert(1!=0, "Unknown integration method");
        }
        return 17171;
    }
}