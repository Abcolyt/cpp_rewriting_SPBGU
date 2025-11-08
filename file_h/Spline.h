#pragma once
#include <array>
#include <algorithm>
#include <utility> 
#include "../file_h/polynomial.h"
#include "../file_h/matrix.h"

uint64_t factorial(const int n)
{
    uint64_t f = 1;
    for (int i = 1; i <= n; ++i)
        f *= i;
    return f;
}

extern enum class output_mode;
template<typename T>
std::vector<T> convert_pairs_to_vector(
    const std::vector<std::pair<T, T>>& input,
    bool take_first_element = true 
) {
    std::vector<T> output;
    output.reserve(input.size());

    
    auto extractor = [take_first_element](const std::pair<T, T>& pair) {
        return take_first_element ? pair.first : pair.second;
    };

    std::transform(
        input.begin(),
        input.end(),
        std::back_inserter(output),
        extractor
    );

    return output;
}

template<typename T>
class Spline {
private:
    std::vector<T> sections;
    std::vector<polynomial<T>> polinoms; // Intervals: sections.size()-1

public:
    
    Spline(T left_border, T right_border, size_t num_sections);
    Spline(std::vector<T> vec);
    ~Spline();

    polynomial<T> getpol(uint64_t index) const{ return polinoms[index]; }
    Spline setpol(uint64_t index, polynomial<T> pol) { polinoms[index] = pol; return *this; }

    uint64_t sections_size()const {
        return sections.size();
    }
    const T& operator[](size_t index) const {
        return sections[index];
    }

    T operator()(const T& x) const {
        for (size_t i = 0; i < sections.size() - 1; ++i) {
            if ( sections[i] <= x  && x <= sections[i + 1]) {
                return polinoms[i](x);
            }
        }
        
        if (x <= sections.front() ) {
            return polinoms.front()(x);
        }
        
        if ( sections.back()<= x) {
            return polinoms[sections.size() - 2](x);
        }
        throw std::out_of_range("x is outside the spline range");
    }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Spline<U>& spline);

    template<uint64_t M_, uint64_t P, typename T>
    friend void Spline_interpolator(const std::vector<std::pair<T, T>>& Array_xy);
};
template<typename T>
Spline<T>::Spline(T left_border, T right_border, size_t num_sections) {
    sections.resize(num_sections);
    polinoms.resize(num_sections-1);
    
    if (left_border >= right_border) {
        throw std::invalid_argument("Left border must be less than right border");
    }
    sections[0] = left_border;
    sections[num_sections - 1] = right_border;
}

template<typename T>
Spline<T>::Spline(std::vector<T> vec) {
    uint64_t Size = vec.size();
    std::sort(vec.begin(), vec.end());
    sections.resize(Size);
    polinoms.resize(Size);
    
    if (vec.size() <=0 ) {
        throw std::invalid_argument("Vector size must match the number of sections");
    }


    for (size_t i = 1; i < Size; ++i) {
        if (vec[i] <= vec[i - 1]) {
            throw std::invalid_argument("Sections must be in ascending order");
        }
    }

    std::copy(vec.begin(), vec.end(), sections.begin());
}

template<typename T>Spline<T>::~Spline()
{
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Spline<T>& spline) {
    os << "Spline with " << spline.sections_size() << " intervals:\n";
    for (size_t i = 0; i < spline.sections_size() -1; ++i) {
        os << "  Interval [" << spline.sections[i]
            << ", " << spline.sections[i + 1] << "): "
            << spline.polinoms[i] << "\n";
    }
    return os;
}

#define SPLINE_LOGS 0
    template<uint64_t M_, uint64_t P_= M_ - 1, typename T>
    Spline<T> Spline_interpolator(std::vector<std::pair<T, T>> Array_xy) {
    #if SPLINE_LOGS==1
        std::cout << "\n";
        for (auto& I : Array_xy) { std::cout << "x:" << I.first << "y:" << I.second << "\n"; }
    #endif
        const uint64_t M = M_+1; 
        const uint64_t P = M_ - 1;
        const size_t N = Array_xy.size();
        const size_t segments = N - 1;
        const size_t coefficients = (M) * segments;

        // Calculation of the total number of equations
        const size_t interpolation_eq = 2 * segments;
        const size_t smoothness_eq = P * (segments - 1);
        const size_t boundary_eq = /*2 **/ P;
        const size_t equations = interpolation_eq + smoothness_eq + boundary_eq;

        if (N == 1) {
            Spline<T> ans(convert_pairs_to_vector(Array_xy));
        
            ans.setpol(0, polynomial<T>(Array_xy[0].second));
            return ans;
        }

        matrix<T> a(equations, coefficients);
        matrix<T> b(equations, 1);

        size_t eq = 0;

        // 1. Interpolation conditions
        for (size_t i = 0; i < segments; ++i) {
            const auto& [x1, y1] = Array_xy[i];
            const auto& [x2, y2] = Array_xy[i + 1];

            // Left point
            for (size_t j = 0; j <M; ++j) {
                a[eq][i * (M) + j] = std::pow(x1, j);
            }
            b[eq][0] = y1;
            eq++;

            // Right point точка
            for (size_t j = 0; j < M; ++j) {
                a[eq][i * M + j] = std::pow(x2, j);
            }
            b[eq][0] = y2;
            eq++;
        }

        // 2. Smoothness conditions (derivatives from 1 to P)
        for (size_t i = 1; i < segments; ++i) {
            const T x = Array_xy[i].first;

            for (size_t p = 1; p <= P; ++p) {
                for (size_t j = p; j < M; ++j) {
                    const T deriv_coeff = factorial(j) / factorial(j - p);
                    const T deriv_value = deriv_coeff * std::pow(x, j - p);

                    a[eq][(i - 1) * M + j] = deriv_value;
                    a[eq][i * M + j] = -deriv_value;
                }
                b[eq][0] = 0;
                eq++;
            }
        }
    #if 0
        // // // //
        // 3. Boundary conditions (derivatives from 1 to P at the ends)
        for (size_t p = 1; p <= P; ++p) {
            // left border
            const T x_start = Array_xy[0].first;
            for (size_t j = p; j < M; ++j) {
                const T deriv_coeff = factorial(j) / factorial(j - p);
                a[eq][0 * M + j] = deriv_coeff * std::pow(x_start, j - p);
            }
           // b[eq][0] = 0;
            eq++;

    #if 0
            // right border
            const double x_end = Array_xy.back().first;
            const size_t last_segment = segments - 1;
            for (size_t j = p; j < M; ++j) {
                const double deriv_coeff = factorial(j) / factorial(j - p);
                a[eq][last_segment * M + j] = deriv_coeff * std::pow(x_end, j - p);
            }
            //b[eq][0] = 0;
            eq++;
    #endif            
        }
        // // // //
    #else
        size_t left_conditions = P / 2;          // number of conditions at the left end
        size_t right_conditions = P - left_conditions; //  right end

        // derivatives from 1 to left_conditions
        for (size_t p = 1; p <= left_conditions; ++p) {
            const T x_start = Array_xy[0].first;
            for (size_t j = p; j < M; ++j) {
                const T deriv_coeff = factorial(j) / factorial(j - p);
                a[eq][0 * M + j] = deriv_coeff * std::pow(x_start, j - p);
            }
            b[eq][0] = 0;  
            eq++;
        }

        // derivatives from 1 to right_conditions
        for (size_t p = 1; p <= right_conditions; ++p) {
            const T x_end = Array_xy.back().first;
            const size_t last_segment = segments - 1;
            for (size_t j = p; j < M; ++j) {
                const T deriv_coeff = factorial(j) / factorial(j - p);
                a[eq][last_segment * M + j] = deriv_coeff * std::pow(x_end, j - p);
            }
            b[eq][0] = 0;  
            eq++;
        }
    #endif
    
    #if SPLINE_LOGS==1
        a.set_output_mode(output_mode::ABBREVIATED);
        std::cout << a << "\n\n";
        std::cout << a.determinant() << "\n\n";
        std::cout<<"b:\n" << b << "\n\n";
    #endif
        matrix<T> x = matrixfunction::solve_system(a, b);
       
    
        Spline<T> ans(convert_pairs_to_vector(Array_xy));
        for (uint64_t i = 0; i < x.getcol() / M; i++) {
            polynomial<T> pol; pol.newsize(M);
            for (uint64_t j = 0; j < M; j++) {
            
                pol[j]=x[M*i+j][0];
            }
            ans.setpol(i,pol);
        }
        return ans;
    }