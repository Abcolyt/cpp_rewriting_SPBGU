﻿#pragma once
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

    //get n normalized newton cotes coefficients(i.e. for [a=0,b=1]
    //get c[j] : Integral[a=0,b=1]F(x)dx of function =~= Sum(j=1,n){c[j]*F(x_j)} 
    template<typename T>
    matrix<T> GetNormalizedСoefficentsNewtonCotes(int n) {
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

            //std::cout << "A:" << A << "\n";
            //std::cout << "b:" << b << "\n";
            matrix<T> coeff = matrixfunction::solve_system(A, b);
            //std::cout << "с_i++:" << A.inverse_M()*b << "\n";
            
            //std::cout << "c_i:" << coeff << '\n';
            return coeff;
        }
    }

#if 0
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
#endif
    namespace special {
        //beta,logbeta,incompletebeta
        namespace betas {
            ////beta function
            // Calculating the logarithm of the beta function: ln(B(p, q))
            long double log_beta_function(long double p, long double q) {
                return std::lgamma(p) + std::lgamma(q) - std::lgamma(p + q);
            }
            //Calculating the beta function
            double beta_function(double p, double q) {
                return std::tgamma(p) * std::tgamma(q) / std::tgamma(p + q);
            }
            // Binomial coefficient calculation
            template<typename T = int>
            T binomial_coefficient(int n, int k) {
                if (k < 0 || k > n) return T(0);
                if (k == 0 || k == n) return T(1);

                T result = 1;
                for (int i = 1; i <= k; ++i) {
                    result *= T(n - i + 1);
                    result /= T(i);
                }
                return result;
            }
            // Forward declaration for incomplete beta function
            template<typename T>
            T incomplete_beta(T z, T a, T b) {
                static_assert(std::is_floating_point_v<T>, "T must be floating point");

                // Boundary cases
                if (z <= T(0)) return T(0);
                if (z >= T(1)) return std::beta(a, b);

                if (a <= T(0) || b <= T(0)) {
                    return std::numeric_limits<T>::quiet_NaN();
                }

                if (z > T(0.5)) {
                    return std::beta(a, b) - incomplete_beta(T(1) - z, b, a);
                }

                const T epsilon = std::numeric_limits<T>::epsilon() * T(10);
                std::cout << "eps[" << epsilon << "]\n";
                const int max_iterations = 100'000;

                T sum = T(0);
                T term = T(1);

                // First term (n=0)
                T current_term = std::pow(z, a) / a;
                sum = current_term;

                for (int n = 1; n < max_iterations; ++n) {
                    // Calculate coefficient using recurrence relation
                    T coeff = (T(n) - b) / T(n) * (a + T(n) - T(1)) / (a + T(n));
                    current_term *= coeff * z;

                    T previous_sum = sum;
                    sum += current_term;

                    // Check 
                    if (std::abs(current_term) < epsilon * std::abs(sum) ||
                        std::abs(sum - previous_sum) < epsilon * std::abs(sum)) {
                        //std::cout << "abs exit\n";
                        break;
                    }
                }

                return sum;
            }
        }
        using namespace betas;

        //IntegralManager() - main function
        namespace euler_type{
            //Calculation of the integral I = ∫[a, b] x ^ j / ((x - a) ^ α * (b - x) ^ β) dx
            long double compute_integral_optimized(long double a, long double b,
                long double alpha, long double beta, int j) {

                if (a >= b) throw std::invalid_argument("a must be less than b");
                if (j < 0) throw std::invalid_argument("j must be non-negative");

                // Special case: a = 0
                if (a == 0.0L) {
                    long double C1 = std::pow(b - a, 1 - alpha - beta + j);
                    return C1 * std::exp(log_beta_function(j - alpha + 1, 1 - beta));
                }

                if (j <= 30) {
                    long double C1 = std::pow(b - a, 1 - alpha - beta);
                    long double C2 = std::pow(a, j);
                    long double C3 = (b - a) / a;

                    long double sum = 0.0L;
                    long double binom_coeff = 1.0L;
                    long double C3_power = 1.0L;

                    for (int k = 0; k <= j; ++k) {
                        long double beta_val = std::exp(log_beta_function(k - alpha + 1, 1 - beta));
                        //std::cout << k - alpha + 1 <<" "<< 1 - beta <<" " << beta_val << '\n';
                        sum += binom_coeff * C3_power * beta_val;

                        if (k < j) {
                            binom_coeff *= static_cast<long double>(j - k) / (k + 1);
                        }
                        C3_power *= C3;
                    }

                    return C1 * C2 * sum;
                }
                else {
                    // The general case: a ≠ 0
                    long double log_C1 = (1 - alpha - beta) * std::log(b - a);
                    long double log_C2 = j * std::log(a);
                    long double log_C3 = std::log(b - a) - std::log(a);

                    std::vector<long double> log_terms(j + 1);

                    long double log_binom = 0.0L; // ln(C(j,0)) = ln(1) = 0

                    for (int k = 0; k <= j; ++k) {
                        // The logarithm of the binomial coefficient + k*ln(C3)
                        long double log_binom_C3 = log_binom + k * log_C3;

                        // log beta function
                        long double log_beta_val = log_beta_function(k - alpha + 1, 1 - beta);

                        // summ: ln(term_k) = ln(binom) + k*ln(C3) + ln(B(...))
                        log_terms[k] = log_binom_C3 + log_beta_val;

                        if (k < j) {
                            // ln(C(j,k+1)) = ln(C(j,k)) + ln(j-k) - ln(k+1)
                            log_binom += std::log(j - k) - std::log(k + 1);
                        }
                    }

                    long double max_log = log_terms[0];
                    for (int k = 1; k <= j; ++k) {
                        if (log_terms[k] > max_log) {
                            max_log = log_terms[k];
                        }
                    }

                    // summarising exp(log_terms[k] - max_log) and multiply on exp(max_log)
                    long double sum_exp = 0.0L;
                    for (int k = 0; k <= j; ++k) {
                        sum_exp += std::exp(log_terms[k] - max_log);
                    }

                    long double log_result = log_C1 + log_C2 + max_log + std::log(sum_exp);
                    return std::exp(log_result);
                }
            }

            ////////
            //Forward declarations for specialized functions
            template<typename TResult, typename TArg>
            TResult Compute_t0_eq_a_T_eq_b_integral(TArg a, TArg b, TArg alpha, TArg beta, int j) { 
                
                if (a >= b) throw std::invalid_argument("a must be less than b");
                if (j < 0) throw std::invalid_argument("j must be non-negative");

                // Special case: a = 0
                if (a == TArg(0)) {
                    TResult C1 = std::pow(b - a, TArg(1) - alpha - beta + j);
                    TResult result = C1 * std::exp(log_beta_function(j - alpha + TArg(1), TArg(1) - beta));
                    return static_cast<TResult>(result);
                }

                if (j <= 30) {
                    TResult C1 = std::pow(b - a, TArg(1) - alpha - beta);
                    TResult C2 = std::pow(a, j);
                    TResult C3 = (b - a) / a;

                    TResult sum = TResult(0);
                    TArg binom_coeff = TArg(1);
                    TArg C3_power = TArg(1);

                    for (int k = 0; k <= j; ++k) {
                        TArg beta_val = std::exp(log_beta_function(
                            k - alpha + TArg(1),
                            TArg(1) - beta
                        ));

                        sum += static_cast<TResult>(binom_coeff * C3_power * beta_val);

                        if (k < j) {
                            binom_coeff *= static_cast<TArg>(j - k) / (k + 1);
                        }
                        C3_power *= C3;
                    }

                    TResult result = C1 * C2 * sum;
                    return result;
                }
                else {
                    // The general case: a ≠ 0
                    TResult log_C1 = (TArg(1) - alpha - beta) * std::log(b - a);
                    TResult log_C2 = j * std::log(a);
                    TResult log_C3 = std::log(b - a) - std::log(a);

                    std::vector<TArg> log_terms(j + 1);
                    TArg log_binom = TArg(0); // ln(C(j,0)) = ln(1) = 0

                    for (int k = 0; k <= j; ++k) {
                        TArg log_binom_C3 = log_binom + k * log_C3;
                        TArg log_beta_val = log_beta_function(
                            k - alpha + TArg(1),
                            TArg(1) - beta
                        );

                        log_terms[k] = log_binom_C3 + log_beta_val;

                        if (k < j) {
                            log_binom += std::log(static_cast<TArg>(j - k)) - std::log(static_cast<TArg>(k + 1));
                        }
                    }

                    TArg max_log = *std::max_element(log_terms.begin(), log_terms.end());

                    TResult sum_exp = TResult(0);
                    for (const auto& term : log_terms) {
                        sum_exp += static_cast<TResult>(std::exp(term - max_log));
                    }

                    TResult log_result = log_C1 + log_C2 + max_log + std::log(sum_exp);
                    TResult result = std::exp(log_result);
                    return result;
                }
                
            }
            //Int[t0;T] x^j dx
            template<typename TResult, typename TArg>
            TResult ComputePolynomialIntegral(TArg t0, TArg T, TArg j) {
                if (j == -1) {
                    // Special case: integral of 1/x = ln(T) - ln(t0)
                    return static_cast<TResult>(std::log(std::abs(T)) - std::log(std::abs(t0)));
                }
                else {
                    auto first = std::pow(T, j + 1);
                    auto second = (t0 != 0 ? std::pow(t0, j + 1) : 0);
                    auto ans = static_cast<TResult>(first - second) / (j + 1);
                    return ans;
                }
            }
            //Int[t0;T] x^j / (b-x^beta))dx
            template<typename TResult, typename TArg>
            TResult ComputeAlphaZeroIntegral(TArg t0, TArg T, TArg b, TArg beta, int j) {
                // Check for valid integration bounds
                if (t0 >= T) {
                    return static_cast<TResult>(0);
                }

                // Check if T > b, which would make the integral invalid
                if (T > b) {
                    throw std::invalid_argument("T must be <= b");
                }

                // Compute shifted variables u = b - x
                TArg u_lower = b - T;
                TArg u_upper = b - t0;

                // Check for convergence at u=0 (which corresponds to x=b)
                if (u_lower == 0) {
                    // Check if we have a non-removable singularity
                    // The integral converges only if k - beta > -1 for all k
                    for (int k = 0; k <= j; k++) {
                        if (k - beta <= -1) {
                            throw std::domain_error("Integral diverges: non-removable singularity at upper bound (x=b)");
                        }
                    }
                }

                TResult result = 0;

                // Optimized computation: precompute b^j and then divide by b at each step
                TArg b_power;
                if (j == 0) {
                    b_power = static_cast<TArg>(1);
                }
                else {
                    b_power = std::pow(b, j);
                }

                // Compute binomial expansion: (b - u)^j = sum_{k=0}^j [C(j,k) * b^(j-k) * (-u)^k]
                for (int k = 0; k <= j; k++) {
                    // Compute binomial coefficient C(j, k)
                    unsigned long long binom = 1;
                    if (k > 0 && k < j) {
                        // Efficient computation of binomial coefficient
                        int k_temp = k;
                        if (k_temp > j - k_temp) k_temp = j - k_temp;
                        for (int i = 1; i <= k_temp; i++) {
                            binom = binom * (j - i + 1) / i;
                        }
                    }
                    else if (k == j) {
                        binom = 1;
                    }

                    // Compute integral of u^(k - beta) from u_lower to u_upper
                    TArg exponent = k - beta;
                    TResult integral_part = ComputePolynomialIntegral<TResult, TArg>(u_lower, u_upper, exponent);

                    // Add to final result with sign (-1)^k
                    TResult term = static_cast<TResult>(binom) * static_cast<TResult>(b_power) * integral_part;
                    if (k % 2 == 1) {
                        term = -term;
                    }
                    result += term;

                    // Update b_power for next iteration: divide by b (unless we're at the last iteration)
                    if (k < j && b != 0) {
                        b_power /= b;
                    }
                }

                return result;
            }
            //Int[t0;T] x^j / ((x-a)^alpha dx
            template<typename TResult, typename TArg>
            TResult ComputeBetaZeroIntegral(TArg t0, TArg T, TArg a, TArg alpha, int j) {
                if (t0 >= T) {
                    return static_cast<TResult>(0);
                }

                TArg u0 = t0 - a;
                TArg U = T - a;

                if (u0 == 0) {
                    // Check if we have a non-removable singularity
                    // The singularity is removable if j >= alpha when alpha is integer,
                    // or more generally, if the limit exists as u->0
                    // For our purposes, we'll check if j >= ceil(alpha) for real alpha
                    if (alpha > j) {
                        throw std::domain_error("Integral diverges: non-removable singularity at lower bound");
                    }
                    // If alpha is integer and j >= alpha, the singularity is removable
                    // If alpha is not integer, we need to be more careful, but for now
                    // we'll assume it's OK if j >= alpha
                }

                TResult result = 0;

                // Compute binomial expansion: (u + a)^j = sum_{k=0}^j [C(j,k) * a^(j-k) * u^k]
                for (int k = 0; k <= j; k++) {
                    // Compute binomial coefficient C(j, k)
                    unsigned long long binom = 1;
                    if (k > 0 && k < j) {
                        // Efficient computation of binomial coefficient
                        int k_temp = k;
                        if (k_temp > j - k_temp) k_temp = j - k_temp;
                        for (int i = 1; i <= k_temp; i++) {
                            binom = binom * (j - i + 1) / i;
                        }
                    }
                    else if (k == j) {
                        binom = 1;
                    }

                    // Compute a^(j-k)
                    TArg a_power;
                    if (j - k == 0) {
                        a_power = static_cast<TArg>(1);
                    }
                    else {
                        a_power = std::pow(a, j - k);
                    }

                    // Compute integral of u^(k - alpha) from u0 to U
                    TArg exponent = (k - alpha);
                    TResult integral_part = ComputePolynomialIntegral<TResult, TArg>(u0, U, exponent);

                    // Add to final result
                    result += static_cast<TResult>(binom) * static_cast<TResult>(a_power) * integral_part;
                }

                return result;
            }
            //Int[t0;T] x^j/((x-a)^alpha * (b-x^beta))dx
            template<typename TResult, typename TArg>
            TResult ComputeGeneralIntegral(TArg t0, TArg T, TArg a, TArg b, TArg alpha, TArg beta, int j) {
                static_assert(std::is_floating_point_v<TArg>, "TArg must be floating point");
                static_assert(std::is_floating_point_v<TResult>, "TResult must be floating point");

                // Input validation
                if (a >= b) {
                    throw std::invalid_argument("a must be less than b");
                }
                if (t0 < a || T > b || t0 > T) {
                    throw std::invalid_argument("Invalid integration limits: a <= t0 <= T <= b required");
                }
                if (alpha >= 1 || beta >= 1) {
                    throw std::invalid_argument("alpha and beta must be less than 1 for convergence");
                }

                // constants
                TArg C1 = std::pow(a, j);
                TArg C2 = std::pow(b - a, TArg(1) - alpha - beta);
                TArg C3 = (b - a) / a;

                TArg t1 = (t0 - a) / (b - a);
                TArg t2 = (T - a) / (b - a);

                TResult sum = TResult(0);

                for (int k = 0; k <= j; ++k) {
                    TArg p = TArg(k) - alpha + TArg(1);
                    TArg q = TArg(1) - beta;

                    TArg binom = binomial_coefficient<TArg>(j, k);

                    TResult beta_diff = static_cast<TResult>(incomplete_beta(t2, p, q) - incomplete_beta(t1, p, q));

                    sum += static_cast<TResult>(binom * std::pow(C3, k)) * beta_diff;
                }

                return static_cast<TResult>(C1 * C2) * sum;
            }
            //////// 
            
            ////
            // Main integrate manager function:  Int[t0;T] x^j/((x-a)^alpha * (b-x)^beta) dx
            template<typename TResult, typename TArg>
            TResult IntegralManager(TArg t0, TArg T, TArg a, TArg b, TArg alpha, TArg beta, int j) {
                // Input validation
                if ((a >= b && beta != 0)) {
                    throw std::invalid_argument("a must be less than b");
                }
                if (t0 < a || (T > b && beta != 0) || t0 > T) {
                    throw std::invalid_argument("Invalid integration limits: a <= t0 <= T <= b required");
                }
                if (!(alpha < 1)) {
                    throw std::invalid_argument("Invalid alpha argument. alpha should be less 1");
                }
                if (!(beta < 1)) {
                    throw std::invalid_argument("Invalid beta argument. beta should be less 1");
                }


                if (t0 == a && T == b) {
                    return Compute_t0_eq_a_T_eq_b_integral<TResult, TArg>(a, b, alpha, beta, j);
                }
                if (a == 0) {
                    return ComputePolynomialIntegral<TResult, TArg>(t0, T, j - alpha);
                }
                if (b == 0) {
                    return ComputePolynomialIntegral<TResult, TArg>(t0, T, j - beta);
                }
                if (alpha == 0 && beta == 0) {
                    return ComputePolynomialIntegral<TResult, TArg>(t0, T, j);
                }

                if (alpha == 0 && beta != 0) {
                    return ComputeAlphaZeroIntegral<TResult, TArg>(t0, T, b, beta, j);
                }

                if (alpha != 0 && beta == 0) {
                    return ComputeBetaZeroIntegral<TResult, TArg>(t0, T, a, alpha, j);
                }

                // Case 5: General case - use full implementation with incomplete beta functions
                return ComputeGeneralIntegral<TResult, TArg>(t0, T, a, b, alpha, beta, j);
            }
            ////
        }
        using namespace euler_type;
        namespace auxiliary {
#define integrate_version 2
#if integrate_version == 1
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
#elif integrate_version==2
        
        template<typename TResult, typename TArg>
        std::unordered_map<int, TResult> fill_integral_map(TArg t0, TArg T,
            TArg a, TArg b,
            TArg alpha, TArg beta,
            int n,
            std::function<TResult(TArg t0, TArg T, TArg a, TArg b, TArg alpha, TArg beta, int)> feature_function) {

            std::unordered_map<int, TResult> integral_map;

            // for j = 0,..., 2*n-1
            for (int j = 0; j < 2 * n; ++j) {
                try {
                    TResult value = feature_function(t0,T,a, b, alpha, beta, j);
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
#endif

////

#define Debug 0
    //main function
        template<typename T = double>
        std::pair<matrix<T>, matrix<T>> GetNGausCoefficientWithP_x_FunctionConstants(const std::unordered_map<int, T>& miu) {

            int n = miu.size() / 2;
            matrix<T> A(n, n), B(n, 1);
            //auto miu = counting_methods_3::fill_integral_map(a, b, alpha, beta, n);
            //input
        //n - число узлов;
#if Debug == 1
            std::cout << "n=" << n << '\n';
            std::cout << "a=" << a << '\n';
            std::cout << "b=" << b << '\n';
            std::cout << "alpha=" << alpha << '\n';
            std::cout << "beta=" << beta << '\n';

#endif



            for (uint64_t i = 0; i < n; i++)
            {
                B[i][0] = -miu.at(n + i);
                for (uint64_t j = 0; j < n; j++)
                {
                    A[i][j] = miu.at(i + j);
                }
            }

#if Debug == 1
            std::cout << "A" << A << '\n';
            std::cout << "B" << B << '\n';
#endif
            matrix<T> coeff = matrixfunction::solve_system(A, B);
#if Debug == 1
            std::cout << "coeff" << coeff << '\n';
            std::cout << "A*coeff" << A * coeff << '\n';
#endif





            polynomial<T> polyn(1);
            for (uint64_t i = 0; i < n; i++)
            {
                //std::cout << "i=" << i << "  (coeff[i][0])" << (coeff[i][0]) << "\n";
                polyn = (polyn >> 1);
                polyn[0] = (coeff[n - 1 - i][0]);


            }
#if Debug == 1
            std::cout << "polyn:" << polyn << '\n';

#endif
            auto roots = polyn.plnm_roots();
            std::vector<T> nodes_in_the_integration_gap;

            int size = roots.size();
            matrix<T>Coeff(size, 1), A2(size), B2(size, 1);

            for (int i = 0; i < size; i++) {
#if Debug == 1

                std::cout << "x_i:" << roots[i].first << '\n';
#endif

                nodes_in_the_integration_gap.push_back(roots[i].first);
                Coeff[i][0] = roots[i].first;
            }
#if Debug == 1
            std::cout << "Coeff" << Coeff << '\n';
#endif  



            //
            //matrix<T>A2(size), B2(size, 1);
            //std::cout << "3:" << nodes_in_the_integration_gap << '\n';
            for (int i = 0; i < size; i++) {
                A2[0][i] = 1;//nodes_in_the_integration_gap[i];
                B2[i][0] = miu.at(i);
                for (int j = 1; j < size; j++) {
                    A2[j][i] = nodes_in_the_integration_gap[i] * A2[j - 1][i];
                }
            }

#if Debug == 1

            std::cout << "A':" << A2 << "\n";
            std::cout << "B':" << B2 << "\n";
#endif
            matrix<T> coeff2 = matrixfunction::solve_system(A2, B2);
            //std::cout << "с'_j:" << A.inverse_M()*b << "\n";
#if Debug == 1

            std::cout << "c'_j:" << coeff2 << '\n';
#endif  


            return { Coeff, coeff2 };
        }

#define Debug 1

        //coeff - nodes, coeff2 - weights
        template<typename TResult, typename TArg, typename T_type = double>
        std::pair<matrix<T_type>, matrix<T_type>> GetNGausCoefficientWithP_x_FunctionConstants(uint64_t n, double t0, double T, double a, double b, double alpha, double beta, std::function<TResult(TArg, TArg, TArg, TArg, TArg, TArg, int)> feature_function = {}) {
            //return GetNGausCoefficientWithP_x_FunctionConstants(counting_methods_3::fill_integral_map(a, b, alpha, beta, n));

#if 1
        //input
    //n - число узлов;
#if Debug == 1
            std::cout << "n=" << n << '\n';
            std::cout << "a=" << a << '\n';
            std::cout << "b=" << b << '\n';
            std::cout << "alpha=" << alpha << '\n';
            std::cout << "beta=" << beta << '\n';

#endif

            matrix<double> A(n, n), B(n, 1);
            auto miu = fill_integral_map(t0, T, a, b, alpha, beta, n, feature_function);


            for (uint64_t i = 0; i < n; i++)
            {
                B[i][0] = -miu[n + i];
                for (uint64_t j = 0; j < n; j++)
                {
                    A[i][j] = miu[i + j];
                }
            }

#if Debug == 1
            std::cout << "A" << A << '\n';
            std::cout << "B" << B << '\n';
#endif
            matrix<double> coeff = matrixfunction::solve_system(A, B);
#if Debug == 1
            std::cout << "coeff" << coeff << '\n';
            std::cout << "A*coeff" << A * coeff << '\n';
#endif





            polynomial<double> polyn(1);
            for (uint64_t i = 0; i < n; i++)
            {
                //std::cout << "i=" << i << "  (coeff[i][0])" << (coeff[i][0]) << "\n";
                polyn = (polyn >> 1);
                polyn[0] = (coeff[n - 1 - i][0]);


            }
#if Debug == 1
            std::cout << "polyn:" << polyn << '\n';

#endif
            auto roots = polyn.plnm_roots();
            std::vector<long double> nodes_in_the_integration_gap;

            int size = roots.size();
            matrix<T_type>Coeff(size, 1), A2(size), B2(size, 1);

            for (int i = 0; i < size; i++) {
#if Debug == 1

                std::cout << "x_i:" << roots[i].first << '\n';
#endif

                nodes_in_the_integration_gap.push_back(roots[i].first);
                Coeff[i][0] = roots[i].first;
            }
#if Debug == 1
            std::cout << "Coeff" << Coeff << '\n';
#endif  



            //
            //matrix<T>A2(size), B2(size, 1);
            //std::cout << "3:" << nodes_in_the_integration_gap << '\n';
            for (int i = 0; i < size; i++) {
                A2[0][i] = 1;//nodes_in_the_integration_gap[i];
                B2[i][0] = miu[i];
                for (int j = 1; j < size; j++) {
                    A2[j][i] = nodes_in_the_integration_gap[i] * A2[j - 1][i];
                }
            }

#if Debug == 1

            std::cout << "A':" << A2 << "\n";
            std::cout << "B':" << B2 << "\n";
#endif
            matrix<T_type> coeff2 = matrixfunction::solve_system(A2, B2);
            //std::cout << "с'_j:" << A.inverse_M()*b << "\n";
#if Debug == 1

            std::cout << "c'_j:" << coeff2 << '\n';
#endif  


            return { Coeff, coeff2 };
#endif
        }


        }
        using namespace auxiliary;
    }
    using namespace special;
    using namespace auxiliary;

#if 0
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
#endif

    

    enum class IntegrateMethod {    
        LEFT_RECTANGLE=0,
        MIDDLE_RECTANGLE,
        TRAPEZOID,
        SIMPSON,
        NEWTON_COTES_3_POINT=4,
        NEWTON_COTES_4_POINT=5,
        NEWTON_COTES_5_POINT,
        NEWTON_COTES_6_POINT,
        NEWTON_COTES_7_POINT,
        NEWTON_COTES_8_POINT,
        NEWTON_COTES_9_POINT=10,
        GAUS_3_POINT,
        GAUS_4_POINT,
        GAUS_5_POINT,
        GAUS_6_POINT
    };



    /*template<typename TResult, typename TArg, IntegrateMethod Method>
    TResult integrate(std::function<TResult(TArg)> f, TArg a, TArg b, int points_of_division, PFeature p_feature ={0,0}) {
        if (p_feature == PFeature{ 0, 0 })return integrate(f,  a,  b, points_of_division);
        if (a > b)throw("a>b");
        TResult sum = 0;
        TResult h = (b - a) / static_cast<TArg>(points_of_division);
        auto a_ = a,b_=a_+h;
        for (uint64_t i = 0; i < points_of_division; i++)
        {
            auto{ dots,veight } = GetNGausCoefficientWithP_x_FunctionConstants(Method, a_, b_, p_feature.alpha, p_feature.beta);
            for (int j = 0; j < dots.getcol(); j++) {
                sum += f(dots[0][j]) * veight[0][j];
            }
            a_ += h, b_ += b_ + h;
        }


    }*/



   //main integrate function. 
    //if (p_feature == { 0, 0 })return integrate(f,  a,  b, int points_of_division);
    //if (a > b)throw("a>b");
    struct PFeature {
        double alpha, beta;
    };
#define integrate_version 2
#if integrate_version == 1
    template<typename TResult, typename TArg, IntegrateMethod Method>
    TResult integrate(std::function<TResult(TArg)> f, TArg a, TArg b, uint64_t interval_of_division, PFeature p_feature = { 0,0 }) {
        TResult sum = 0;
        if (a > b){ throw("a>b"); }

        TArg h = (b - a) / static_cast<TArg>(interval_of_division);

        if constexpr (Method == IntegrateMethod::LEFT_RECTANGLE) {
            if (p_feature.alpha == 0)sum += f(a);
            for (int i = 1; i < interval_of_division; ++i) {
                sum += f(a + i * h);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::MIDDLE_RECTANGLE) {
            for (int i = 1; i < interval_of_division; ++i) {
                sum += f(a + i * h - h / 2);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::TRAPEZOID) {
            if (p_feature.alpha == 0) { sum += f(a)/2; }
            if (p_feature.beta == 0) { sum += f(b)/2; }
            for (int i = 1; i < interval_of_division; ++i) {
                sum += f(a + i * h);
            }
            return sum * h;
        }
        else if constexpr (Method == IntegrateMethod::SIMPSON) {
            if(p_feature.alpha == 0){ sum+=f(a); }
            if(p_feature.beta  == 0){ sum+=f(b); }
            TResult sum_x4 = 0;
            TResult sum_x2 = 0;
            for (int i = 1; i < interval_of_division; i += 2) {
                sum_x4 += f(a + i * h);
            }
            for (int i = 2; i < interval_of_division; i += 2) {
                sum_x2 += f(a + i * h);
            }
            sum += 4 * sum_x4 + 2 * sum_x2;
            return sum * h / 3;
        }
        else if constexpr ((static_cast<int>(IntegrateMethod::NEWTON_COTES_3_POINT) <= static_cast<int>(Method)) && (static_cast<int>(Method) <= static_cast<int>(IntegrateMethod::NEWTON_COTES_9_POINT))) {
            int n = static_cast<int>(Method);
            TResult micro_h = h / static_cast<TArg>(n-1);
            matrix<TResult> M = GetNormalizedСoefficentsNewtonCotes<TResult>(n);
            if (p_feature.alpha== 0)sum += M[0][0] * f(a);
            for (int j = 0; j < interval_of_division; j++)
            {
                for (int i = 1; i < n-1; ++i) {

                    sum += M[0][i] * f(a + h * j + micro_h * i);
                }
            }
            for (int j = 1; j < interval_of_division; j++)
            {
                sum += 2*M[0][0] * f(a + h * j);
            }
            if (p_feature.beta == 0)sum += M[0][0] * f(b);

            return sum * h;
        }
        else if constexpr ((static_cast<int>(IntegrateMethod::GAUS_3_POINT) <= static_cast<int>(Method)) && (static_cast<int>(Method) <= static_cast<int>(IntegrateMethod::GAUS_6_POINT))) {
            auto a_ = a, b_ = a_ + h;
            auto miu_coeff = counting_methods_3::fill_integral_map(0, h, p_feature.alpha, p_feature.beta, static_cast<int>(Method) - 8);
            //std::cout << a;
            for (uint64_t i = 0; i < interval_of_division; i++)
            {
#if 1
                std::pair<matrix<TArg>, matrix<TResult>> M = GetNGausCoefficientWithP_x_FunctionConstants(static_cast<int>(Method) - 8,a_, b_, p_feature.alpha, p_feature.beta);
#else
                std::pair<matrix<TArg>, matrix<TResult>> M = GetNGausCoefficientWithP_x_FunctionConstants(miu_coeff);
#endif
                auto dots = M.first, weight = M.second;
                for (int j = 0; j < dots.getcol(); j++) {
                    sum += f(dots[0][j]) * weight[0][j];
                }
                a_ += h, b_ += h;
            }
            return sum;
        }

        else {
            static_assert(1!=0, "Unknown integration method");
        }
        return 404;
    }
#else

template<typename TResult, typename TArg, IntegrateMethod Method>
TResult integrate(std::function<TResult(TArg)> f,
    TArg a, TArg b, uint64_t interval_of_division, PFeature p_feature = { 0,0 },
    std::function<TResult(TArg t0, TArg T, TArg a, TArg b, TArg alpha, TArg beta, int)> feature_function = euler_type::IntegralManager<TResult, TArg>) {
    TResult sum = 0;
    bool sh_be_reversed = 0;
    if (a > b) {
        throw("a>b");
        std::swap(a, b); 
        sh_be_reversed = 1;
    }
    TArg h = (b - a) / static_cast<TArg>(interval_of_division);
    if constexpr (Method == IntegrateMethod::LEFT_RECTANGLE) {
        if (p_feature.alpha == 0)sum += f(a);
        for (int i = 1; i < interval_of_division; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    }
    else if constexpr (Method == IntegrateMethod::MIDDLE_RECTANGLE) {
        for (int i = 1; i < interval_of_division; ++i) {
            sum += f(a + i * h - h / 2);
        }
        return sum * h;
    }
    else if constexpr (Method == IntegrateMethod::TRAPEZOID) {
        if (p_feature.alpha == 0) { sum += f(a) / 2; }
        if (p_feature.beta == 0) { sum += f(b) / 2; }
        for (int i = 1; i < interval_of_division; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    }
    else if constexpr (Method == IntegrateMethod::SIMPSON) {
        if (p_feature.alpha == 0) { sum += f(a); }
        if (p_feature.beta == 0) { sum += f(b); }
        TResult sum_x4 = 0;
        TResult sum_x2 = 0;
        for (int i = 1; i < interval_of_division; i += 2) {
            sum_x4 += f(a + i * h);
        }
        for (int i = 2; i < interval_of_division; i += 2) {
            sum_x2 += f(a + i * h);
        }
        sum += 4 * sum_x4 + 2 * sum_x2;
        return sum * h / 3;
    }
    else if constexpr ((static_cast<int>(IntegrateMethod::NEWTON_COTES_3_POINT) <= static_cast<int>(Method)) && (static_cast<int>(Method) <= static_cast<int>(IntegrateMethod::NEWTON_COTES_9_POINT))) {
        int n = static_cast<int>(Method);
        TResult micro_h = h / static_cast<TArg>(n - 1);
        matrix<TResult> M = GetNormalizedСoefficentsNewtonCotes<TResult>(n);
        if (p_feature.alpha == 0)sum += M[0][0] * f(a);
        for (int j = 0; j < interval_of_division; j++)
        {
            for (int i = 1; i < n - 1; ++i) {

                sum += M[0][i] * f(a + h * j + micro_h * i);
            }
        }
        for (int j = 1; j < interval_of_division; j++)
        {
            sum += 2 * M[0][0] * f(a + h * j);
        }
        if (p_feature.beta == 0)sum += M[0][0] * f(b);

        return sum * h;
    }
    else if constexpr ((static_cast<int>(IntegrateMethod::GAUS_3_POINT) <= static_cast<int>(Method)) && (static_cast<int>(Method) <= static_cast<int>(IntegrateMethod::GAUS_6_POINT))) {
        auto a_ = a, b_ = a_ + h;
        //std::cout << a;
        for (uint64_t i = 0; i < interval_of_division; i++)
        {
#if 0
            std::pair<matrix<TArg>, matrix<TResult>> M = GetNGausCoefficientWithP_x_FunctionConstants(static_cast<int>(Method) - 8, a_ ,b_,a,b, p_feature.alpha, p_feature.beta, feature_function);
#else
            auto miu_coeff = counting_methods_3::fill_integral_map(a_,b_,a,b, p_feature.alpha, p_feature.beta, static_cast<int>(Method) - 8, feature_function);
            std::pair<matrix<TArg>, matrix<TResult>> M = GetNGausCoefficientWithP_x_FunctionConstants(miu_coeff);
#endif
            auto dots = M.first, weight = M.second;
            for (int j = 0; j < dots.getcol(); j++) {
                sum += f(dots[0][j]) * weight[0][j];
            }
            a_ += h, b_ += h;
        }
        return sum;
    }

    else {
        static_assert(1 != 0, "Unknown integration method");
    }
    return 404;
}

#endif

}