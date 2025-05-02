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
    std::vector<polynomial<T>> polinoms; // Интервалов: sections.size()-1

public:
    
    Spline(T left_border, T right_border, size_t num_sections);
    Spline(const std::vector<T>& vec);
    ~Spline();
    uint64_t sections_size()const {
        return sections.size();
    }
    const T& operator[](size_t index) const {
        return sections[index];
    }

    T operator()(const T& x) const {
        for (size_t i = 0; i < sections.size() - 1; ++i) {
            if (x >= sections[i] && x < sections[i + 1]) {
                return polinoms[i](x);
            }
        }
        
        if (x == sections.back()) {
            return polinoms.back()(x);
        }
        throw std::out_of_range("x is outside the spline range");
    }


    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Spline<U>& spline);
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
Spline<T>::Spline(const std::vector<T>& vec) {
    uint64_t Size = vec.size();

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

template<uint64_t M, uint64_t P, typename T>void Spline_build(const std::vector<std::pair<T,T>> Array_xy) {
    //using N = Array_xy.size();
    Spline<T> spline(convert_pairs_to_vector(Array_xy));
    matrix<double> a;
    a.set_output_mode(output_mode::ABBREVIATED);
    a.zeros(M*Array_xy.size(), M*Array_xy.size());
    for(uint64_t i0 = 0; i0 < M*Array_xy.size(); i0+=2)
    {
        for (uint64_t i = 0; i < Array_xy.size(); i++)
        {
            if (i0 * M + i < M*Array_xy.size()) {
                a[i0][i0*M+i] = factorial() ;
            }
            
        }
    }

    std::cout << a;
    matrix<double> b(Array_xy.size(), 1);
    for(uint64_t i = 0; i < Array_xy.size(); i++)
    {
        b[i][0] = Array_xy[i].second;
    }
    std::cout << b;
}
