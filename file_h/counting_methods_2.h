#pragma once
#include <vector>
#include <algorithm>
#include "file_h/polynomial.h"
#include <random>
#include <functional>
#include <memory>

#include <cmath>
#include <corecrt_math_defines.h>
#include <iomanip> 

namespace counting_methods_2 {

    namespace Polynomial_interpolation {
        namespace nuton {
            //
            /* polynomial<double> W(k,) */
            template<typename P>class Nuton :public polynomial<P>
            {
            public:
                Nuton(const std::vector<std::pair<P, P>>& Array_xy);
                Nuton();
                ~Nuton();
                P nuton_interpolation(uint64_t k);
            private:
                std::vector<std::pair<P, P>> Array_xy;

                std::vector<P> divided_difference_order_0_to_n;
                //razdelenay raznost
                

                polynomial<P> W_k(uint64_t k, P x);
            };
            template<typename P>polynomial<P> Nuton<P>::W_k(uint64_t k, P x) {
            
            
            }
            template<typename P>Nuton<P>::Nuton(const std::vector<std::pair<P, P>>& Array_xy) {
                this->Array_xy = Array_xy;
                divided_difference_order_0_to_n.reserve(Array_xy.size());
            }
            template<typename P>Nuton<P>::Nuton() : polynomial<P>()
            {
            }
            template<typename P>Nuton<P>::~Nuton()
            {
            }

            template<typename P>P Nuton<P>::nuton_interpolation(uint64_t k) {
                std::transform(
                    Array_xy.begin(), Array_xy.end(), // out
                    std::back_inserter(divided_difference_order_0_to_n), // in
                    [](const std::pair<P, P>& pair) { return pair.second; } // f()
                );
                #define LOGS  for (uint64_t i = 0; i < Array_xy.size(); i++) \
                {                                                         \
                    std::cout << divided_difference_order_0_to_n[i];     \
                }   

                LOGS;
                for (size_t i = 1; i < Array_xy.size(); i++)
                {   
                    divided_difference_order_0_to_n[i - 1] = (divided_difference_order_0_to_n[i] - divided_difference_order_0_to_n[i - 1]) / ((Array_xy[i]).first - Array_xy[i - 1].first);
                }
                std::cout << "\n\n\n\n\n\n\n\n";
                


                LOGS;
                return 0;
            }
        
        }

        namespace nuton2 {

            template<typename P>std::vector<std::pair<P, P>> extractUniqueY(std::vector<std::pair<P, P>> Array_xy) {

                std::sort(
                    Array_xy.begin(),
                    Array_xy.end(),
                    [](const auto& a, const auto& b) {return a.first < b.first; }
                );

                std::vector<std::pair<P, P>> result;
                result.push_back(Array_xy[0]);
                for (size_t i = 1; i < Array_xy.size(); ++i) {
                    // Пропускаем пары, где x совпадает с предыдущим
                    if (Array_xy[i].first == Array_xy[i - 1].first) {

                        continue;
                    }
                    result.push_back(Array_xy[i]);
                }
                return result;
            }

            template<typename P>polynomial<P> nuton_interpolation(std::vector<std::pair<P, P>> Array_xy) {

                Array_xy = extractUniqueY(Array_xy);

                //#define LOGS  for (uint64_t i = 0; i < Array_xy.size(); i++) \
                //{                                                         \
                //    std::cout << Array_xy[i].first<<"  "<< Array_xy[i].second <<"\n";     \
                //} 

                //LOGS;

                polynomial<P> Ans = 0, w_k = 1;

                for (size_t i = 1; i < Array_xy.size(); i++)
                {
                    /*if ((Array_xy[0].second) == 0)
                    {
                        return Ans;
                    }*/

                    // std::cout << "Ans:" << Ans << " (Array_xy[0].second)=( " << (Array_xy[0].second) << ")*w_k " << w_k << '\n';
                    Ans = Ans + w_k * (Array_xy[0].second);
                    for (size_t j = 0; j < Array_xy.size() - i; j++)
                    {
                        Array_xy[j].second = (Array_xy[j + 1].second - Array_xy[j].second) / (Array_xy[j + i].first - Array_xy[j].first);
                    }


                    //std::cout << "Ans:" << Ans << "\n";
                    w_k = w_k * ((static_cast<polynomial<P>>(1) >> 1) - Array_xy[i - 1].first);
                    //std::cout << "i"<< i<< '\n';
                    //LOGS
                }
                Ans = Ans.cutbag();
                return Ans;
            }

            polynomial<int> generateRandomIntCoefficients(int min_degree = 3, int max_degree = 10, double min_c = -10, double max_c = 10) {
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


            //generate with polinomial function
            template<typename T>std::vector<std::pair<T, T>> generatePointsFuncPtr(int k, T x0, T step, polynomial<T>& polynom, const std::function<T(polynomial<T>, T)>& F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < k; ++i) {
                    T x = x0 + i * step;
                    points.emplace_back(x, F(polynom, x));
                }
                return points;
            }
            template<typename T, typename Func>std::vector<std::pair<T, T>> generatePoints_equally_sufficient(int k, T x0, T step, Func F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < k; ++i) {
                    T x = x0 + i * step;
                    points.emplace_back(x, F(x));
                }
                return points;
            }

            template<typename T, typename Func>std::vector<std::pair<T, T>> generatePoints_optimal(int n, T a, T b, Func F) {
                std::vector<std::pair<T, T>> points;
                for (int i = 0; i < n; ++i) {
                    T x = (0.5) * ((b - a) * std::cos(M_PI * ((2.0 * i + 1) / (2 * (n)))) + (b + a));
                    //std::cout << "\n" << x<<" F(x):"<<F(x);
                    points.emplace_back(x, F(x));
                }
                return points;
            }

            //4-working function
            template<typename P, typename Func>polynomial<P> N_n(int n, P a, P b, Func F) {
                return nuton_interpolation(extractUniqueY(generatePoints_equally_sufficient(n, a, (b - a) / n, F)));
            }
            //4-working function
            template<typename P, typename Func>polynomial<P> N_optn(int n, P a, P b, Func F) {
                return nuton_interpolation(extractUniqueY(generatePoints_optimal(n, a, b, F)));
            }

            //F := working function
            template<typename P, typename Func>double RN_n(int n, P a, P b, Func F) {
                auto table = extractUniqueY(generatePoints_equally_sufficient(n, a, (b - a) / n, F));
                auto interpolinom = nuton_interpolation(table);
                double Max_ans = 0;
                for (const auto& p : table) {
                    Max_ans = std::max(Max_ans, std::abs(polynomialfunctions::f_polyn_x0_(interpolinom, p.first) - p.second));
                }
                return Max_ans;
            }

            template<typename P, typename Func>double RN_optn(int n, P a, P b, Func F) {
                auto table = extractUniqueY(generatePoints_optimal(n, a, b, F));
                auto interpolinom = nuton_interpolation(table);
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


            //Lagranggg
            template<typename P>polynomial<P> w_k0(std::vector<P> array, int T0)
            {
            polynomial<P> a(1);
            for (uint64_t i = 0; i < array.size(); i++)
            {
                if (i != T0) {
                    polynomial<P> b(1);
                    b = 1; b = b >> 1; 
                    //std::cout << b << '\n';
                    b = b - array[i];
                    //std::cout << b << '\n';
                    a = a * b;
                    //std::cout << a << '\n';
                }
            }
            return a;
            }

            template<typename P, typename Func>polynomial<P> L_n(int n, P a, P b, Func F);

        }




        namespace Lagrang {

            //prohodit chereze zadannie tochki
            //djbavit zavisimost ot naklona s bokov
            void Lagrang() {}
        }
    }
    namespace Spline_interpolation {}
}
