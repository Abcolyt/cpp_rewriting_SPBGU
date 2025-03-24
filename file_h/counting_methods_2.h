#pragma once
#include <vector>
#include <algorithm>
#include "file_h/polynomial.h"

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
                P divided_difference0_to_k(uint64_t k);
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

            template<typename P>P Nuton<P>::divided_difference0_to_k(uint64_t k) {
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

            template<typename P>
            std::vector<std::pair<P, P>> extractUniqueY(std::vector<std::pair<P, P>> Array_xy) {

                std::sort(
                    Array_xy.begin(),
                    Array_xy.end(),
                    [](const auto& a, const auto& b) {
                        return a.first < b.first;
                    }
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

            template<typename P>polynomial<P> divided_difference0_to_k(std::vector<std::pair<P, P>> Array_xy) {
                
                Array_xy = extractUniqueY(Array_xy);

                //#define LOGS  for (uint64_t i = 0; i < Array_xy.size(); i++) \
                //{                                                         \
                //    std::cout << Array_xy[i].first<<"  "<< Array_xy[i].second <<"\n";     \
                //} 

                //LOGS;

                polynomial<P> Ans=0, w_k = 1;

                for (size_t i = 1; i < Array_xy.size(); i++)
                {
                    /*if ((Array_xy[0].second) == 0)
                    {
                        return Ans;
                    }*/

                   // std::cout << "Ans:" << Ans << " (Array_xy[0].second)=( " << (Array_xy[0].second) << ")*w_k " << w_k << '\n';
                    Ans = Ans + w_k * (Array_xy[0].second);
                    for (size_t j = 0; j < Array_xy.size()-i; j++)
                    {
                        Array_xy[j].second = (Array_xy[j + 1].second - Array_xy[j].second) / (Array_xy[j + i].first - Array_xy[j].first);
                    }
                    
                    
                    //std::cout << "Ans:" << Ans << "\n";
                    w_k = w_k * ((static_cast<polynomial<P>>(1) >> 1) - Array_xy[i - 1].first) ;
                    //std::cout << "i"<< i<< '\n';
                    //LOGS
                }
                Ans=Ans.cutbag();
                return Ans;
            }
        }

        namespace Lagrang {

            //prohodit chereze zadannie tochki
            //djbavit zavisimost ot naklona s bokov
            void Lagrang() {}
        }
    }
    namespace Spline_interpolation {}
}
