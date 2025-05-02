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
//rezerved
/*
template<uint64_t M, uint64_t P, typename T>void Spline_build(const std::vector<std::pair<T,T>> Array_xy) {
    //using N = Array_xy.size();
    Spline<T> spline(convert_pairs_to_vector(Array_xy));
    matrix<double> a;
    a.set_output_mode(output_mode::ABBREVIATED);
    a.zeros(M*Array_xy.size(), M*Array_xy.size());
    for(uint64_t i0 = 0; i0 < M*Array_xy.size(); i0+=M)
    {

        for (uint64_t i = 0; i < M; i++)
        {
            if (i0 * M + i < M*Array_xy.size()) {
                a[i0][i0*M+i] = factorial(i)*std::pow(Array_xy[i0].first,i);
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

*/

//rezerved_2
/*
template<uint64_t M, uint64_t P, typename T>
void Spline_build(const std::vector<std::pair<T, T>>& Array_xy) {
    const size_t N = Array_xy.size();
    const size_t total_coeffs = M * (N - 1); // Коэффициенты для каждого интервала
    const size_t equations = M * (N - 1);     // Общее количество уравнений

    matrix<double> a;
    a.set_output_mode(output_mode::ABBREVIATED);
    a.zeros(equations, equations);

    std::vector<double> b(equations, 0.0);

    // Условия интерполяции (2 уравнения на каждый интервал)
    size_t eq = 0;
    for (size_t i = 0; i < N - 1; ++i) {
        double x_left = Array_xy[i].first;
        double x_right = Array_xy[i + 1].first;

        // Левая точка интервала
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x_left, j);
        }
        b[eq] = Array_xy[i].second;
        eq++;

        // Правая точка интервала
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x_right, j);
        }
        b[eq] = Array_xy[i + 1].second;
        eq++;
    }

    // Условия гладкости (производные до порядка P)
    for (size_t i = 1; i < N - 1; ++i) {
        double x_common = Array_xy[i].first;

        for (size_t p = 1; p <= P; ++p) {
            // Производная от левого полинома
            for (size_t j = p; j < M; ++j) {
                a[eq][(i - 1) * M + j] = factorial(j) / factorial(j - p) * std::pow(x_common, j - p);
            }

            // Производная от правого полинома (с противоположным знаком)
            for (size_t j = p; j < M; ++j) {
                a[eq][i * M + j] = -(factorial(j) / factorial(j - p) * std::pow(x_common, j - p));
            }

            eq++;
        }
    }

    // Граничные условия (пример: естественный сплайн)
    if (eq < equations) {
        // Вторая производная на левом конце равна 0
        size_t i = 0;
        for (size_t j = 2; j < M; ++j) {
            a[eq][i * M + j] = factorial(j) / factorial(j - 2) * std::pow(Array_xy[i].first, j - 2);
        }
        eq++;

        // Вторая производная на правом конце равна 0
        i = N - 2;
        for (size_t j = 2; j < M; ++j) {
            a[eq][i * M + j] = factorial(j) / factorial(j - 2) * std::pow(Array_xy.back().first, j - 2);
        }
        eq++;
    }




    // Решение системы уравнений
    matrix<double> x = a.inverse_M()*b;
    std::cout << x;
    // Далее используем x для установки коэффициентов сплайна
}
*/

//main rezerved
/*
template<uint64_t M, uint64_t P, typename T>void Spline_build(const std::vector<std::pair<T, T>>& Array_xy) {
    const size_t N = Array_xy.size();
    const size_t segments = N - 1;
    const size_t coefficients = M * segments;
    const size_t interpolation_eq = 2 * segments;
    const size_t smoothness_eq = P * (segments - 1);
    const size_t boundary_eq = 2; // Пример: естественный сплайн
    const size_t equations = interpolation_eq + smoothness_eq + boundary_eq;

    matrix<double> a(equations, equations);
    matrix<double> b(equations, 1); // Правильный размер: M*(N-1) x 1

    // Заполнение матрицы a и вектора b
    size_t eq = 0;

    // 1. Условия интерполяции (2*(N-1) уравнений)
    for (size_t i = 0; i < N - 1; ++i) {
        double x_left = Array_xy[i].first;
        double x_right = Array_xy[i + 1].first;

        // Левая точка сегмента
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x_left, j);
        }
        b[eq][0] = Array_xy[i].second;
        eq++;

        // Правая точка сегмента
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x_right, j);
        }
        b[eq][0] = Array_xy[i + 1].second;
        eq++;
    }

    // 2. Условия гладкости (P*(N-2) уравнений)
    for (size_t i = 1; i < N - 1; ++i) {
        double x = Array_xy[i].first;

        for (size_t p = 1; p <= P; ++p) {
            // Производные левого и правого полиномов
            for (size_t j = p; j < M; ++j) {
                double deriv_coeff = factorial(j) / factorial(j - p);
                a[eq][(i - 1) * M + j] = deriv_coeff * std::pow(x, j - p);
                a[eq][i * M + j] = -deriv_coeff * std::pow(x, j - p);
            }
            b[eq][0] = 0; // Условие: производные совпадают
            eq++;
        }
    }

    // 3. Граничные условия (оставшиеся уравнения)
    // Пример: естественный сплайн (вторая производная на концах = 0)
    for (size_t p = 2; p <= 2; ++p) { // Вторая производная
        // Левый конец
        double x = Array_xy[0].first;
        for (size_t j = p; j < M; ++j) {
            a[eq][0 * M + j] = factorial(j) / factorial(j - p) * std::pow(x, j - p);
        }
        b[eq][0] = 0;
        eq++;

        // Правый конец
        x = Array_xy.back().first;
        size_t last_segment = N - 2;
        for (size_t j = p; j < M; ++j) {
            a[eq][last_segment * M + j] = factorial(j) / factorial(j - p) * std::pow(x, j - p);
        }
        b[eq][0] = 0;
        eq++;
    }
    a.set_output_mode(output_mode::ABBREVIATED);
    std::cout << a<<"\nzer good\n";
    // Решение системы
    matrix<double> x = a.inverse_M() * b; // Если матрица квадратная и невырожденная
    // ИЛИ использовать метод solve(a, b) с проверкой на вырожденность
    std::cout << x;
}
*/

//rezerv_ XZ
/*
template<uint64_t M, uint64_t P, typename T>
void Spline_build(const std::vector<std::pair<T, T>>& Array_xy) {
    const size_t N = Array_xy.size();
    const size_t segments = N - 1;
    const size_t coefficients = M * segments;
    const size_t equations = 2 * segments + P * (segments - 1) + 2;

    matrix<double> a(equations, coefficients);
    matrix<double> b(equations, 1);

    size_t eq = 0;

    // 1. Интерполяционные условия
    for (size_t i = 0; i < segments; ++i) {
        const auto& [x1, y1] = Array_xy[i];
        const auto& [x2, y2] = Array_xy[i + 1];

        // Левая точка
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x1, j);
        }
        b[eq][0] = y1;
        eq++;

        // Правая точка
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x2, j);
        }
        b[eq][0] = y2;
        eq++;
    }

    // 2. Условия гладкости (первые производные)
    for (size_t i = 1; i < segments; ++i) {
        const double x = Array_xy[i].first;

        for (size_t j = 1; j < M; ++j) {
            const double deriv = j * std::pow(x, j - 1);
            a[eq][(i - 1) * M + j] = deriv;
            a[eq][i * M + j] = -deriv;
        }
        b[eq][0] = 0;
        eq++;
    }

    // 3. Граничные условия (вторая производная = 0)
    const double x_start = Array_xy[0].first;
    for (size_t j = 2; j < M; ++j) {
        a[eq][0 * M + j] = j * (j - 1) * std::pow(x_start, j - 2);
    }
    b[eq][0] = 0;
    eq++;

    const double x_end = Array_xy.back().first;
    for (size_t j = 2; j < M; ++j) {
        a[eq][(segments - 1) * M + j] = j * (j - 1) * std::pow(x_end, j - 2);
    }
    b[eq][0] = 0;
    eq++;
    a.set_output_mode(output_mode::ABBREVIATED);
    std::cout << a << "\nzer good\n";

    // Решение системы
    matrix<double> x = a.inverse_M() * b;
    std::cout << "Coefficients:\n" << x;
}
*/

template<uint64_t M_, uint64_t P, typename T>
void Spline_build(const std::vector<std::pair<T, T>>& Array_xy) {
    const uint64_t M = M_+1; 

    const size_t N = Array_xy.size();
    const size_t segments = N - 1;
    const size_t coefficients = (M) * segments;

    // Расчет общего количества уравнений
    const size_t interpolation_eq = 2 * segments;
    const size_t smoothness_eq = P * (segments - 1);
    const size_t boundary_eq = 2 * P;
    const size_t equations = interpolation_eq + smoothness_eq + boundary_eq;

    matrix<double> a(equations, coefficients);
    matrix<double> b(equations, 1);

    size_t eq = 0;

    // 1. Условия интерполяции
    for (size_t i = 0; i < segments; ++i) {
        const auto& [x1, y1] = Array_xy[i];
        const auto& [x2, y2] = Array_xy[i + 1];

        // Левая точка
        for (size_t j = 0; j <M; ++j) {
            a[eq][i * (M) + j] = std::pow(x1, j);
        }
        b[eq][0] = y1;
        eq++;

        // Правая точка
        for (size_t j = 0; j < M; ++j) {
            a[eq][i * M + j] = std::pow(x2, j);
        }
        b[eq][0] = y2;
        eq++;
    }

    // 2. Условия гладкости (производные от 1 до P)
    for (size_t i = 1; i < segments; ++i) {
        const double x = Array_xy[i].first;

        for (size_t p = 1; p <= P; ++p) {
            for (size_t j = p; j < M; ++j) {
                const double deriv_coeff = factorial(j) / factorial(j - p);
                const double deriv_value = deriv_coeff * std::pow(x, j - p);

                a[eq][(i - 1) * M + j] = deriv_value;
                a[eq][i * M + j] = -deriv_value;
            }
            b[eq][0] = 0;
            eq++;
        }
    }

    // 3. Граничные условия (производные от 1 до P на концах)
    if (boundary_eq > 0) {
        // Левая граница (вторая производная = 0)
        const double x_start = Array_xy[0].first;
        for (size_t j = 2; j < M; ++j) {
            a[eq][0 * M + j] = j * (j - 1) * std::pow(x_start, j - 2);
        }
        b[eq][0] = 0;
        eq++;

        // Правая граница (вторая производная = 0)
        const double x_end = Array_xy.back().first;
        const size_t last_segment = segments - 1;
        for (size_t j = 2; j < M; ++j) {
            a[eq][last_segment * M + j] = j * (j - 1) * std::pow(x_end, j - 2);
        }
        b[eq][0] = 0;
        eq++;
    }
    a.set_output_mode(output_mode::ABBREVIATED);
    std::cout << a << "\nzer good\n";
    //std::cout << a.determinant() << "\nzer good\n";
    //// Решение системы
    //matrix<double> x = a.inverse_M() * b;
    //std::cout << "Coefficients:\n" << x;
}