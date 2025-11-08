#pragma once
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <memory>
#include <iomanip> 
#include <unordered_map>
#include <map>
#include <cmath>
#include <corecrt_math_defines.h>

#include <numeric>

#include "core/polynomial.h"
#include "core/Array_xy_To_.h"
#include "interpolation/Spline.h"

namespace counting_methods_2 {

    namespace Polynomial_interpolation {
      

        namespace nuton2 {

            enum class methodical_error
            {
                nuton_interpolation_double_generatePoints_equally_sufficient_ = 0,
                nuton_interpolation_double_generatePoints_optimal = 1,
                Lagrang_interpolation_double_generatePoints_equally_sufficient_ = 2,
                Lagrang_interpolation_double_generatePoints_optimal = 3,
                Alternativ_Lagrang_interpolation_double_generatePoints_equally_sufficient_=4,
                Alternativ_Lagrang_interpolation_double_generatePoints_optimal=5,
                Alternativ_nuton_interpolation_double_generatePoints_equally_sufficient_ = 6,
                Alternativ_nuton_interpolation_double_generatePoints_optimal = 7,
                


            };
            using enum methodical_error;
            const std::map<std::string, enum class methodical_error> teoretical_max_error{
                {"nuton_interpolation<double>generatePoints_equally_sufficient_", nuton_interpolation_double_generatePoints_equally_sufficient_},
                {"nuton_interpolation<double>generatePoints_optimal", nuton_interpolation_double_generatePoints_optimal},

                {"Lagrang_interpolation<double>generatePoints_equally_sufficient_",Lagrang_interpolation_double_generatePoints_equally_sufficient_},
                {"Lagrang_interpolation<double>generatePoints_optimal",Lagrang_interpolation_double_generatePoints_optimal},

                {"Alternativ_Lagrang_interpolation<double>generatePoints_equally_sufficient_",Alternativ_Lagrang_interpolation_double_generatePoints_equally_sufficient_},
                { "Alternativ_Lagrang_interpolation<double>generatePoints_optimal",Alternativ_Lagrang_interpolation_double_generatePoints_optimal },

                {"Alternativ_nuton_interpolation<double>generatePoints_equally_sufficient_",Alternativ_nuton_interpolation_double_generatePoints_equally_sufficient_},
                { "Alternativ_nuton_interpolation<double>generatePoints_optimal",Alternativ_nuton_interpolation_double_generatePoints_equally_sufficient_ }
            };


            uint64_t factorial(uint64_t n) { 
                return n>=1?factorial(n-1)*n:1; 
            }

            template<class T>
            inline constexpr T pow(const T base, unsigned const exponent)
            {
                return (exponent == 0) ? 1 : (base * pow(base, exponent - 1));
            }

            template<typename P, typename FuncInterpolation>P methodic_error(enum methodical_error type, const std::vector<std::pair<P, P>>& Array_xy,FuncInterpolation func_interpolation) {
                using enum methodical_error;
                
                polynomial<P> pol= func_interpolation(Array_xy);
                std::vector <P> array_x;
                
                std::transform(Array_xy.begin(), Array_xy.end(),
                    std::back_inserter(array_x),
                    [](const auto& pair) { return pair.first; });
                pol.the_root_constructor(array_x);
                /*for (auto& i : Array_xy){ array_x.push_back(i.first); }*/
                           
                auto min = *std::min_element(array_x.begin(), array_x.end());
                //std::cout << " " << (pol) <<"  "<< (min)<<"  (" << pol(min) << " )";
                auto max = *std::max_element(array_x.begin(), array_x.end());
                //std::cout << (int)type << "\n";
                switch(type)
                {
                case Alternativ_Lagrang_interpolation_double_generatePoints_equally_sufficient_:
                    return (pol.maximum_abs(min, max) / factorial(Array_xy.size() + 1));
                    break;
                case Alternativ_Lagrang_interpolation_double_generatePoints_optimal:
                    return 2 * (std::pow((max - min) / 4, Array_xy.size() + 1)) / factorial(Array_xy.size() + 1) ;
                    break;

                case Lagrang_interpolation_double_generatePoints_equally_sufficient_:      
                    return (pol.maximum_abs(min, max) / factorial(Array_xy.size() + 1));
                    break;
                case Lagrang_interpolation_double_generatePoints_optimal:
                    return 2 * (std::pow((max - min) / 4, Array_xy.size() + 1)) / factorial(Array_xy.size() + 1);
                    break;

                case nuton_interpolation_double_generatePoints_equally_sufficient_:        break;
                case nuton_interpolation_double_generatePoints_optimal:                    break;
                default:return P(1);
                }
                return P(1);
            }

            template<typename P>std::vector<std::pair<P, P>> filter_by_unique_x(std::vector<std::pair<P, P>> Array_xy) {

                std::sort(
                    Array_xy.begin(),
                    Array_xy.end(),
                    [](const auto& a, const auto& b) {return a.first < b.first; }
                );

                std::vector<std::pair<P, P>> result;
                result.push_back(Array_xy[0]);
                for (size_t i = 1; i < Array_xy.size(); ++i) {
                   
                    if (Array_xy[i].first == Array_xy[i - 1].first) {

                        continue;
                    }
                    result.push_back(Array_xy[i]);
                }
                return result;
            }

            polynomial<int> generateRandomIntCoefficients(int min_degree = 3, int max_degree = 10, double min_c = -10, double max_c = 10)
            {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dist_deg(min_degree, max_degree);
                polynomial<int> Ans;
                Ans.newsize(dist_deg(gen));

                std::uniform_int_distribution<> dist_c(min_c, max_c);
                for (int i = 0; i < Ans.get_deg(); ++i) {
                    Ans[i] = (dist_c(gen));
                }

                Ans.cutbag();

                return Ans;
            }

            template<typename P>polynomial<P> nuton_interpolation(std::vector<std::pair<P, P>> Array_xy) {

                Array_xy = filter_by_unique_x(Array_xy);
                polynomial<P> Ans = 0, w_k = 1;

                for (size_t i = 1; i < Array_xy.size(); i++)
                {

                    Ans = Ans + w_k * (Array_xy[0].second);
                    for (size_t j = 0; j < Array_xy.size() - i; j++)
                    {
                        Array_xy[j].second = (Array_xy[j + 1].second - Array_xy[j].second) / (Array_xy[j + i].first - Array_xy[j].first);
                    }

                    w_k = w_k * ((static_cast<polynomial<P>>(1) >> 1) - Array_xy[i - 1].first);
                }
                Ans = Ans.cutbag();
                return Ans;
            }

            template<typename P>polynomial<P> Alternativ_nuton_interpolation(std::vector<std::pair<P, P>> Array_xy) {
                Array_xy = filter_by_unique_x(Array_xy);
                polynomial<P> Ans = 0, w_k = 1;
                size_t n = Array_xy.size();

                for (size_t i = 0; i < n; ++i) {
                    
                    P diff = 0;
                    for (size_t j = 0; j <= i; ++j) {
                        
                        P denominator = 1;
                        for (size_t k = 0; k <= i; ++k) {
                            if (k != j) {
                                denominator *= (Array_xy[j].first - Array_xy[k].first);
                            }
                        }
                        diff += Array_xy[j].second / denominator;
                    }

                    
                    Ans = Ans + w_k * diff;

                    
                    if (i < n - 1) {
                        w_k = w_k * ((static_cast<polynomial<P>>(1) >> 1) - Array_xy[i].first);
                    }
                }

                return Ans.cutbag();
            }

            //generate with points for interpolinomial function
            //==========
            template<typename T>std::vector<std::pair<T, T>> generatePointsFuncPtr(int k, T x0, T step, polynomial<T>& polynom, const std::function<T(polynomial<T>, T)>& F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < k; ++i) {
                    T x = x0 + i * step;
                    points.emplace_back(x, F(polynom, x));
                }
                return points;
            }
            template<typename T, typename Func>std::vector<std::pair<T, T>> generatePoints_equally_sufficient_with_step_size(int k, T x0, T step, Func F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < k; ++i) {
                    T x = x0 + i * step;
                    points.emplace_back(x, F(x));
                }
                return points;
            }

            template<typename T, typename Func>std::vector<std::pair<T, T>> generatePoints_equally_sufficient_(int k, T a, T b, Func F) {
                T step = (b - a)/k;
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < k; ++i) {
                    T x = a + i * step;
                    points.emplace_back(x, F(x));
                }
                return points;
            }

            template<typename T, typename Func>std::vector<std::pair<T, T>> generatePoints_optimal(int n, T a, T b, Func F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < n; ++i) {
                    T x = (0.5) * ((b - a) * std::cos(M_PI * ((2.0 * i + 1) / (2 * (n)+1))) + (b + a));
                    //std::cout << "\n" << x<<" F(x):"<<F(x);
                    points.emplace_back(x, F(x));
                }
                return points;
            }
            //==========
            //4-working function
            template<typename P, typename Func>polynomial<P> N_n(int n, P a, P b, Func F) {
                return nuton_interpolation(filter_by_unique_x(generatePoints_equally_sufficient_with_step_size(n, a, (b - a) / n, F)));
            }
            //4-working function
            template<typename P, typename Func>polynomial<P> N_optn(int n, P a, P b, Func F) {
                return nuton_interpolation(filter_by_unique_x(generatePoints_optimal(n, a, b, F)));
            }

            //F := working function
            template<typename P, typename Func>double RN_n(int n, P a, P b, Func F) {
                auto table = filter_by_unique_x(generatePoints_equally_sufficient_with_step_size(n, a, (b - a) / n, F));
                auto interpolinom = nuton_interpolation(table);
                double Max_ans = 0;
                for (const auto& p : table) {
                    Max_ans = std::max(Max_ans, std::abs(polynomialfunctions::f_polyn_x0_(interpolinom, p.first) - p.second));
                }
                return Max_ans;
            }

            template<typename P, typename Func>double RN_optn(int n, P a, P b, Func F) {
                auto table = filter_by_unique_x(generatePoints_optimal(n, a, b, F));
                auto interpolinom = nuton_interpolation(table);
                double Max_ans = 0;
                for (const auto& p : table) {
                    Max_ans = std::max(Max_ans, std::abs(polynomialfunctions::f_polyn_x0_(interpolinom, p.first) - p.second));
                }
                return Max_ans;
            }


            //Lagranggg
            template<typename P>polynomial<P> w_k_T0(std::vector<P> array, int64_t T0)
            {
            polynomial<P> a(1);
            for (uint64_t i = 0; i < array.size(); i++)
            {
                if (i != T0) {
                    polynomial<P> b(1);
                    b = 1; b = b >> 1; 
                    b = b - array[i];
                    a = a * b;
                }
            }
            return a;
            }
            template<typename P>polynomial<P> w_k_T0(std::vector<std::pair<P,P>>& array_xy, int64_t T0)
            {
                polynomial<P> a(1);
                for (uint64_t i = 0; i < array_xy.size(); i++)
                {
                    if (i != T0) {
                        polynomial<P> b(1);
                        b = 1; b = b >> 1;
                        b = b - array_xy[i].first;
                        a = a * b;
                    }
                }
                return a;
            }

            template<typename P>polynomial<P> Lagrang_interpolation(std::vector<std::pair<P, P>> Array_xy) {
                Array_xy = filter_by_unique_x(Array_xy);
                polynomial<P> ans(0);
                for (uint64_t i = 0; i < Array_xy.size(); i++)
                {
                    auto pl = w_k_T0(Array_xy, i);
                    ans = ans + (pl / (pl(Array_xy[i].first))) * Array_xy[i].second;
                }
                return ans;
            }
            
            template<typename P>polynomial<P> Alternativ_Lagrang_interpolation(std::vector<std::pair<P, P>> Array_xy) {
                Array_xy = filter_by_unique_x(Array_xy);
                
                matrix<P> A(Array_xy.size()),B(Array_xy.size(),1);
                
                for (uint64_t i = 0; i < Array_xy.size(); i++)
                {
                    P el = Array_xy[i].first, el_= 1;
                    for (size_t j = 0; j < Array_xy.size(); j++)
                    {
                        A[i][j] = el_ ;
                        el_ = el_ * el;
                    }
                    B[i][0] = Array_xy[i].second;
                }
                A = (A.inverse_M()) * B;
                //std::cout << A << '\n';
                polynomial<P> ans(0);
                ans.set_deg(A.getcol()) ;
                for (uint64_t i = 0; i < A.getcol()-1; i++)
                {
                    ans[i] = A[i][0];
                    //std::cout << A[i][0] << '\n';
                }
                return ans;
                
            }
#if 0
            template<typename P, typename Func>polynomial<P> L_n(int n, P a, P b, Func F) {
                return Lagrang_interpolation(filter_by_unique_x(generatePoints_equally_sufficient_with_step_size(n, a, (b - a) / n, F)));
            
            }
            template<typename P, typename Func>polynomial<P> L_optn(int n, P a, P b, Func F) {
                return Lagrang_interpolation(filter_by_unique_x(generatePoints_optimal(n, a, b, F)));
            }

            template<typename P, typename Func>double RL_n(int n, P a, P b, Func F) {
                auto table = filter_by_unique_x(generatePoints_equally_sufficient_with_step_size(n, a, (b - a) / n, F));
                auto interpolinom = Lagrang_interpolation(table);
                double Max_ans = 0;
                for (const auto& p : table) {
                    Max_ans = std::max(Max_ans, std::abs(polynomialfunctions::f_polyn_x0_(interpolinom, p.first) - p.second));
                }
                return Max_ans;
            }
            template<typename P, typename Func>double RL_optn(int n, P a, P b, Func F) {
                auto table = filter_by_unique_x(generatePoints_optimal(n, a, b, F));
                auto interpolinom = Lagrang_interpolation(table);
                double Max_ans = 0;
                for (const auto& p : table) {
                    Max_ans = std::max(Max_ans, std::abs(polynomialfunctions::f_polyn_x0_(interpolinom, p.first) - p.second));
                }
                return Max_ans;
            }

            //max n 
            void show_nuton(int n = 10, int m_ = 10) {
                double a = -2 * M_PI, b = 2 * M_PI;

                std::cout << std::left
                    << std::setw(20) << "(n)"
                    << std::setw(30) << "(m)"
                    << std::setw(25) << "RN_n"
                    << std::setw(25) << "RNopt_n"
                    << "\n------------------------------------------------------------\n";

                for (size_t i = 1; i < n + m_; i++)
                {

                    auto F = [](double x) { return std::cos(x) / std::sin(x) + x * x; };
                    auto poly = N_optn(i, a, b, F);

                    // std::cout << "\nPOLINOM N_optn:" << poly;

                    double rn = RN_n(i, a, b, F);
                    double rn_opt = RN_optn(i, a, b, F);

                    std::cout << std::left
                        << std::setw(20) << i
                        << std::setw(30) << i + m_
                        << std::setw(25) << std::setprecision(6) << rn
                        << std::setw(25) << std::setprecision(6) << rn_opt
                        << "\n";
                }




            }
#endif


            //n - degree of interpolation polinom,  m - degree of check 
                template<typename P,
                typename TheTypeOfFunctionThatBeingInterpolated,
                typename TypeOfInterpolationFunc, 
                typename TypeFuncOfPointGeneration = decltype(generatePoints_optimal<P, TheTypeOfFunctionThatBeingInterpolated>) >
                auto error_of_the_interpolation_function(uint64_t n, uint64_t m,P a,P b,
                    TheTypeOfFunctionThatBeingInterpolated F,
                    TypeOfInterpolationFunc interpolation,
                    TypeFuncOfPointGeneration point_generation) -> std::pair<double, decltype(interpolation(std::declval<std::vector<std::pair<P, P>>>() )) >
                {
                    auto table = filter_by_unique_x(point_generation(n+m, a, b, F));

                    auto interpolinom = interpolation(
                        //std::vector<std::pair<P, P>>(table.begin(), table.begin() + std::min(m, table.size()))
                        filter_by_unique_x(point_generation(n, a, b, F))
                    );

                    double Max_ans = 0;
                    for (const auto& p : table) {
                        Max_ans = std::max(
                            Max_ans,
                            std::abs((interpolinom(p.first)) - p.second)
                        );
                    }
                    return { Max_ans,std::move(interpolinom) };
                }
                bool should_add_point(uint64_t i, uint64_t n, uint64_t num_display) {
                    if (num_display == 0) return true;

                    if (n == 0) return false;
                    if (n == 1) return (i == 1);  

                    double step = (num_display > 1) ? (n - 1.0) / (num_display - 1) : 0;

                    if (i == 1 || i == n) return true;

                    if (num_display == 1) return (i == n);

                    double position = (i - 1.0) / step;
                    double fractional = position - std::floor(position);

                    return std::min(fractional, 1.0 - fractional) < 0.01;
                }
            template<typename P=double,typename TypeOfInterpolationFunc /*= decltype(Lagrang_interpolation<double>)*/, typename TheTypeOfFunctionThatBeingInterpolated>
            std::stringstream show_interpolation_statistic(TypeOfInterpolationFunc func_interpolation,int n = 10, int m_ = 99, const std::string& interpolation_name="Unknown",
                TheTypeOfFunctionThatBeingInterpolated F = [](double x) { return std::cos(x) / std::sin(x) + x * x; },
                P a = LDBL_EPSILON -2 * M_PI,P b = 2 * M_PI -LDBL_EPSILON,
                const uint64_t number_of_displaying_function=0
            ) {
                using InterpolationResult = decltype(func_interpolation(std::declval<std::vector<std::pair<P, P>>>()));
                std::stringstream Ans;
                int s1 = std::max(std::toupper(log10(n)), 3), s2 = std::max(std::toupper(log10(n + m_)), 3), s3 = 2 + interpolation_name.size() + 2, s4 = std::string{"R_teoretical_n"  }.size()+6,s5 = std::string{ "R_teoretical_opt__n" }.size()+6;
                Ans << "\n"
                    << std::left
                    << "+" + (std::string(s1, '-') + "+") + (std::string(s2, '-') + "+") + (std::string(s3, '-') + "+") + (std::string(s3 + 5, '-') + "+") + (std::string(s4, '-') + "+") + (std::string(s5, '-'))
                    << "+" << "\n"
                    
                    << std::left
                    << std::setw(1) << "|"
                    << std::setw(s1) << "(n)"
                    << std::setw(1)<<"|"
                    << std::setw(s2) << "(m)"
                    << std::setw(1) << "|"
                    << std::setw(s3) << "R_" + interpolation_name+ "_n"
                    << std::setw(1) << "|"
                    << std::setw(s3+5) << "R_" + interpolation_name + "_opt__n"
                   

                    << std::setw(1) << "|"
                    << std::setw(s4) << "R_teoretical_n/M_n+1"
                    << std::setw(1) << "|"
                    << std::setw(s5) << "R_teoretical_opt__n/M_n+1"
                    << std::setw(1) << "|"

                    << "\n+" + (std::string(s1, '-') + "+") + (std::string(s2, '-') + "+") + (std::string(s3, '-') + "+") + (std::string(s3 + 5, '-') + "+") + (std::string(s4, '-') + "+") + (std::string(s5, '-'))
                    << "+" << "\n";

                std::vector<std::function<P(P)>> functions_n, functions_opt, functions_n_abs_diff, functions_opt_abs_diff;
                //==
                functions_n.reserve(n + 1); 
                functions_n.push_back(F);

                functions_opt.reserve(n + 1); 
                functions_opt.push_back(F);
                //==
                functions_n_abs_diff.reserve(n );
                

                functions_opt_abs_diff.reserve(n );
                
                //==
                for (uint64_t i = 2; i <= n ; i++)
                {
                    auto rn = error_of_the_interpolation_function(i, m_, a, b, F, func_interpolation, generatePoints_equally_sufficient_);
                    auto rn_opt = error_of_the_interpolation_function(i,m_, a, b, F, func_interpolation, generatePoints_optimal);
                    //
                    auto interp_fn = rn.second;//func_interpolation(generatePoints_equally_sufficient_(i, a, b, F));
                    auto interp_fopt = rn_opt.second;//func_interpolation(generatePoints_optimal(i, a, b, F));

                    if (should_add_point(i, n, number_of_displaying_function)) {
                    //auto poly = N_optn(i, a, b, F);
                    functions_n.push_back(rn.second);
                    functions_opt.push_back(rn_opt.second);
                    functions_n_abs_diff.push_back([interp_fn, F](P x) {return std::abs(interp_fn(x) - F(x));});
                    functions_opt_abs_diff.push_back([interp_fopt, F](P x) { return std::abs(interp_fopt(x) - F(x)); });
                    }
                    //

                    Ans << std::left
                        << std::setw(1) << "|"
                        << std::setw(s1) << i
                        << std::setw(1) << "|"
                        << std::setw(s2) << m_
                        << std::setw(1) << "|"
                        << std::setw(s3) << std::setprecision(12) << rn.first
                        << std::setw(1) << "|"
                        << std::setw(s3 + 5) << std::setprecision(12) << rn_opt.first;
                    if constexpr (std::is_same_v<InterpolationResult, polynomial<P>>) {
                        Ans << std::setw(1) << "|"
                            << std::setw(s4) << std::setprecision(12) << methodic_error(teoretical_max_error.at(interpolation_name + "generatePoints_equally_sufficient_"), generatePoints_equally_sufficient_(i, a, b, F), func_interpolation)
                            << std::setw(1) << "|"
                            << std::setw(s5) << std::setprecision(12) << methodic_error(teoretical_max_error.at(interpolation_name + "generatePoints_optimal"), generatePoints_optimal(i, a, b, F), func_interpolation)
                            << std::setw(1) << "|"
                            << "\n";
                    }
                    else if constexpr (std::is_same_v<InterpolationResult, Spline<P>>) {
                        Ans << std::setw(1) << "|"
                            << std::setw(s4) << std::setprecision(12) <<"has not"
                            << std::setw(1) << "|"
                            << std::setw(s5) << std::setprecision(12) << "has not"
                            << std::setw(1) << "|"
                            << "\n";
                    }

                }
                Ans<< "+" + (std::string(s1, '-') + "+") + (std::string(s2, '-') + "+") + (std::string(s3, '-') + "+") + (std::string(s3 + 5, '-') + "+") + (std::string(s4, '-') + "+") + (std::string(s5, '-')) << "+"
                    << "\n";
                std::cout << Ans.str();

                // // // //
                    std::vector<std::function<P(P)>> addition_functions_n, addition_functions_opt, addition_functions_n_abs_diff, addition_functions_opt_abs_diff;
                if (interpolation_name == "Spline_interpolator<3, 2, double>") {


                    auto func_interpolation_two = nuton_interpolation<double>;

                    auto rn = error_of_the_interpolation_function(n, m_, a, b, F, func_interpolation_two, generatePoints_equally_sufficient_);
                    auto rn_opt = error_of_the_interpolation_function(n, m_, a, b, F, func_interpolation_two, generatePoints_optimal);
                    //
                    auto interp_fn = rn.second;//func_interpolation(generatePoints_equally_sufficient_(i, a, b, F));
                    auto interp_fopt = rn_opt.second;//func_interpolation(generatePoints_optimal(i, a, b, F));

                   /* std::cout << rn.second << "\n";
                    std::cout << rn_opt.second << "\n";*/


                    addition_functions_n.push_back(rn.second);
                    addition_functions_opt.push_back(rn_opt.second);
                    addition_functions_n_abs_diff.push_back([interp_fn, F](P x) {return std::abs(interp_fn(x) - F(x));});
                    addition_functions_opt_abs_diff.push_back([interp_fopt, F](P x) { return std::abs(interp_fopt(x) - F(x)); });
                }
                // // // //
#if __has_include(<SFML/Graphics.hpp>)

                draw_functions(functions_n, a, b, std::vector<std::pair<P, P>>{}, ("R_" + interpolation_name + "_n"), addition_functions_n);
                draw_functions(functions_opt, a, b, std::vector<std::pair<P, P>>{}, ("R_" + interpolation_name + "_opt__n"), addition_functions_opt);

                draw_functions(functions_n_abs_diff, a, b, std::vector<std::pair<P, P>>{}, ("R_" + interpolation_name + "_n_abs_diff"), addition_functions_n_abs_diff);
                draw_functions(functions_opt_abs_diff, a, b, std::vector<std::pair<P, P>>{}, ("R_" + interpolation_name + "opt_n_abs_diff"), addition_functions_opt_abs_diff);
                

#else
                std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
                <<  "not working draw_functions(functions_n)" <<"\n"
                << "not working draw_functions(functions_opt)" << "\n";
#endif
                return Ans;
            }
            
}
            };

namespace aproximate {
    // dot of ortagoality
    //degree of polinom
    template<typename P>std::vector<polynomial<P>> orthogonal_polynomials(const std::vector<P>& x, int n) {
        std::vector<polynomial<P>> polys;
        if (n < 0) return polys;
        polys.push_back(polynomial<P>(1));
        if (n == 0) return polys;

        //  alpha1
        P alpha1 = std::accumulate(x.begin(), x.end(), P(0)) / static_cast<P>(x.size());
        
        //P alpha1 = 1;

        //  q1(x) = x - alpha1
        polynomial<P> q1;
        q1.newsize(2); 
        q1[0] = -alpha1;
        q1[1] = 1;
        polys.push_back(q1);


        if (n == 1) return polys;

        for (int j = 1; j < n; j++) {
            std::vector<P> qj_vals, qjm1_vals;
            for (auto xi : x) {
                qj_vals.push_back(polys[j](xi));
                qjm1_vals.push_back(polys[j - 1](xi));
            }

            P num_alpha = 0, den_alpha = 0;
            for (int i = 0; i < x.size(); i++) {
                num_alpha += x[i] * qj_vals[i] * qj_vals[i];
                den_alpha += qj_vals[i] * qj_vals[i];
            }
            P alpha_j1 = num_alpha / den_alpha;

            P num_beta = 0, den_beta = 0;
            for (int i = 0; i < x.size(); i++) {
                num_beta += x[i] * qj_vals[i] * qjm1_vals[i];
                den_beta += qjm1_vals[i] * qjm1_vals[i];
            }
            P beta_j = num_beta / den_beta;

            polynomial<P> x_minus_alpha;
            x_minus_alpha.newsize(2);
            x_minus_alpha[0] = -alpha_j1;
            x_minus_alpha[1] = 1;

            polynomial<P> term1 = x_minus_alpha * polys[j];
            polynomial<P> term2 = polys[j - 1] * beta_j;

            polys.push_back(term1 - term2);
        }

        return polys;
    }
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

    // // // ///
    template<typename P>polynomial<P> least_squares_normal(std::vector<std::pair<P, P>> Array_xy, uint64_t degree_of_the_polynomial = 5) {
        auto V = matrix<double>::vander(convert_pairs_to_vector(Array_xy), degree_of_the_polynomial);
        matrix<double> b((Array_xy.size()), 1);
        int i = 0;
        for (auto& I : convert_pairs_to_vector(Array_xy, false)) {
            b[i][0] = I;
            i++;
        }
        auto matrx_ans = matrixfunction::sanitize_zeros(V.pseudo_inverse() * b, 1e-9);
        polynomial<double> pol_ans; pol_ans.set_deg(degree_of_the_polynomial);
        for (size_t i = 0; i < degree_of_the_polynomial; i++)
        {
            pol_ans[i] = matrx_ans[i][0];
        }
        return pol_ans;
    }

    template<typename P>
    polynomial<P> least_squares_orthogonal(const std::vector<std::pair<P, P>>& points, uint64_t n) {
        std::vector<P> x, y;
        for (const auto& point : points) {
            x.push_back(point.first);
            y.push_back(point.second);
        }

        std::vector<polynomial<P>> phi = counting_methods_2::aproximate::orthogonal_polynomials(x, n);

        std::vector<P> c(n, 0.0);

        for (uint64_t k = 0; k < n; k++) {
            P numerator = 0.0;
            P denominator = 0.0;
            //Numerator: projection of the vector y onto the base vector φₖ.
            //Denominator : the square of the norm φₖ at discrete points.

            for (size_t i = 0; i < x.size(); i++) {
                P phi_k_x = phi[k](x[i]);
                numerator += y[i] * phi_k_x;
                denominator += phi_k_x * phi_k_x;
            }
            //k-th basic vector
            c[k] = numerator / denominator;
        }

        polynomial<P> result;
        for (uint64_t k = 0; k < n; k++) {
            result = result + (phi[k] * c[k]);
        }

        return result;
    }


    // // // ///

    template<typename P>
    P sum_squared_errors(const polynomial<P>& poly,
        const std::vector<std::pair<P, P>>& points) {
        P total_error = 0;

        for (const auto& point : points) {
            P x = point.first;
            //P y_true = ;
            P y_pred = poly(x); 
            P error = y_pred - point.second;
            total_error += error * error;
        }

        return total_error;
    }

    template<typename P = double, typename TheTypeOfFunctionThatBeingInterpolated>
    std::stringstream show_aproximate_statistic(
        TheTypeOfFunctionThatBeingInterpolated F,  
        int n = 10,
        int size = 99,
        P a = LDBL_EPSILON - 2 * M_PI,
        P b = 2 * M_PI - LDBL_EPSILON
    )
    {
        

        std::stringstream Ans;
        const std::string aproximate_name = "least_squares";

        int s1 = std::max(static_cast<int>(std::ceil(std::log10(n))), 3);
        int s3 = 2 + aproximate_name.size() + 2;

        std::string horizontal_line =
            "+" + std::string(s1, '-') + "+" +
            std::string(s3, '-') + "+" +
            std::string(s3 + 5, '-') + "+";

        Ans << "\n"
            << std::left
            << horizontal_line << "\n"
            << std::setw(1) << "|"
            << std::setw(s1) << "(n)"
            << std::setw(1) << "|"
            << std::setw(s3) << "R_" + aproximate_name + "_n"
            << std::setw(1) << "|"
            << std::setw(s3 + 5) << "R_" + aproximate_name + "_opt__n"
            << std::setw(1) << "|\n"
            << horizontal_line << "\n";

        // lines of function
        std::vector<std::function<P(P)>> functions_normal, functions_ortog;
        functions_normal.reserve(n + 1);
        functions_normal.push_back(F);
        functions_ortog.reserve(n + 1);
        functions_ortog.push_back(F);

        // // // //data
        std::vector<std::pair<double, double>> clean_points = counting_methods_2::Polynomial_interpolation::nuton2::generatePoints_equally_sufficient_(size, a, b, F);
        auto noisy_points = addNoiseToPoints(clean_points, 3, 0.2);
        noisy_points.insert(noisy_points.end(), clean_points.begin(), clean_points.end());
        // // // //

        for (uint64_t i = 1; i <= n; i++)
        {

            functions_normal.push_back(least_squares_normal(noisy_points, i));
            functions_ortog.push_back(least_squares_orthogonal(noisy_points, i));
            //
            auto sum_of_squared_errors_norm = sum_squared_errors(least_squares_normal(noisy_points, i), noisy_points);
            auto sum_of_squared_errors_ortog = sum_squared_errors(least_squares_orthogonal(noisy_points, i), noisy_points);

   
            Ans << std::left
                << std::setw(1) << "|"
                << std::setw(s1) << i
                << std::setw(1) << "|"
                << std::setw(s3) << std::setprecision(12) << sum_of_squared_errors_norm
                << std::setw(1) << "|"
                << std::setw(s3 + 5) << std::setprecision(12) << sum_of_squared_errors_ortog
                << std::setw(1) << "|\n";
        }
        Ans << horizontal_line << "\n";
        std::cout << Ans.str();

#if __has_include(<SFML/Graphics.hpp>)
        draw_functions(functions_normal, a, b, noisy_points, aproximate_name + "_normal");
        draw_functions(functions_ortog, a, b, noisy_points,  aproximate_name + "_orthogonal");
#else
        std::cout << "\n__has_include(<SFML/Graphics.hpp>)==0\n"
            << "not working draw_functions(functions_normal)\n"
            << "not working draw_functions(functions_ortog)\n";
#endif
        return Ans;
    }



}

}
