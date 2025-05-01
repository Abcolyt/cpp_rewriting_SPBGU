#pragma once
#include <array>
#include <algorithm>
#include <utility> 
#include "../file_h/polynomial.h"
#include "../file_h/matrix.h"

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

template<typename T, size_t Size>
class Spline {
private:
    std::array<T, Size> sections;
    std::array<polynomial<T>, Size - 1> polinoms; // Интервалов: Size-1

public:

    // Определите конструктор и деструктор
    Spline(T left_border, T right_border);
    Spline(const std::vector<T>& vec);
    ~Spline();

    const T& operator[](size_t index) const {
        return sections[index];
    }

    T operator()(const T& x) const {
        for (size_t i = 0; i < Size - 1; ++i) {
            if (x >= sections[i] && x < sections[i + 1]) {
                return polinoms[i](x);
            }
        }
        
        if (x == sections.back()) {
            return polinoms.back()(x);
        }
        throw std::out_of_range("x is outside the spline range");
    }


    template<typename U, size_t N>
    friend std::ostream& operator<<(std::ostream& os, const Spline<U, N>& spline);
};
template<typename T, size_t Size>
Spline<T,Size>::Spline(T left_border, T right_border) {
    static_assert(Size >= 2, "Spline requires at least 2 sections"); // Проверка на этапе компиляции

    if (left_border >= right_border) {
        throw std::invalid_argument("Left border must be less than right border");
    }
    sections[0] = left_border;
    sections[Size - 1] = right_border;
}

template<typename T, size_t Size>
Spline<T, Size>::Spline(const std::vector<T>& vec) {
    static_assert(Size >= 1, "Spline requires at least 1 section");

    if (vec.size() != Size) {
        throw std::invalid_argument("Vector size must match the number of sections");
    }


    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] <= vec[i - 1]) {
            throw std::invalid_argument("Sections must be in ascending order");
        }
    }

    std::copy(vec.begin(), vec.end(), sections.begin());
}

template<typename T, size_t Size>
Spline<T, Size>::~Spline()
{
}
template<uint64_t M, uint64_t P, typename T>
void Spline_build(const std::vector<std::pair<T,T>> Array_xy) {
    //using N = Array_xy.size();
    Spline<T, Array_xy.size()> spline(convert_pairs_to_vector(Array_xy));
    matrix<double> a;
    a.zeros(Array_xy.size(), Array_xy.size());
    for (uint64_t i = 0; i < 0; i++)
    {

    }

    std::cout << a;
    
}

template<typename T, size_t Size>
std::ostream& operator<<(std::ostream& os, const Spline<T, Size>& spline) {
    os << "Spline with " << Size << " sections:\n";
    for (size_t i = 0; i < Size - 1; ++i) {
        os << "  Interval [" << spline.sections[i]
            << ", " << spline.sections[i + 1] << "): "
            << spline.polinoms[i] << "\n";
    }
    return os;
}