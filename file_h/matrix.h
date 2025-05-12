#pragma once
#include<iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <random>

#include <stdexcept>
#include <cmath>
#include <numeric>
#include <algorithm>


extern enum class output_mode;
template<typename T> class matrix;
template<typename T>struct LUResult {
    matrix<T> L;
    matrix<T> U;
    matrix<T> P; 
};
namespace matrixfunction {
    template<typename T>T power_method(const matrix<T>& A, const matrix<T>& Vec0, double epsilon = 1e-6, int max_iter = 1000);

    template<typename T>T power_method(const matrix<T>& A, double epsilon = 1e-6, int max_iter = 1000) {
        return matrixfunction::power_method(A, matrix<T>::ones(A.getcol(), 1), epsilon, max_iter);
    }

    template <typename T>bool is_equal(const matrix<T>& a, const matrix<T>& b, T epsilon = static_cast<T>(1e-6));

    template <typename T>matrix<int> elementwise_equal_matrix(const matrix<T>& a, const matrix<T>& b, T epsilon = static_cast<T>(1e-6));

#if 1
    // Вспомогательные функции для прямого и обратного хода
    template<typename T>
    matrix<T> forward_substitution(const matrix<T>& L, const matrix<T>& b) {
        int n = L.getrow();
        matrix<T> y(n, 1);
        for (int i = 0; i < n; ++i) {
            T sum = b[i][0];
            for (int j = 0; j < i; ++j) {
                sum -= L[i][j] * y[j][0];
            }
            // Предполагаем, что L имеет единицы на диагонали (как в LU-разложении)
            y[i][0] = sum;
        }
        return y;
    }

    template<typename T>
    matrix<T> backward_substitution(const matrix<T>& U, const matrix<T>& y) {
        int n = U.getrow();
        matrix<T> x(n, 1);
        for (int i = n - 1; i >= 0; --i) {
            T sum = y[i][0];
            for (int j = i + 1; j < n; ++j) {
                sum -= U[i][j] * x[j][0];
            }
            x[i][0] = sum / U[i][i];
        }
        return x;
    }

    // 2. Решение системы с обработкой вырожденных матриц
    template<typename T>
    matrix<T> solve_system(const matrix<T>& A, const matrix<T>& b) {
        auto lup = A.LUP();
        matrix<T> Pb = lup.P * b;
        matrix<T> y = forward_substitution(lup.L, Pb);
        matrix<T> x = backward_substitution(lup.U, y);
        return x;
    }

    // 3. Модифицированный обратный степенной метод
    template<typename T>
    std::pair<T, matrix<T>> inverse_power_method(
        const matrix<T>& A,
        T sigma,
        const matrix<T>& initial_vec,
        double delta = 1e-8,  // Уменьшенный порог для учета малых компонент
        double epsilon = 1e-6,
        int max_iter = 500)   // Увеличенное число итераций
    {
        int n = A.getrow();
        matrix<T> y = initial_vec;
        T norm = y.norm();
        if (norm < 1e-10) throw std::runtime_error("Initial vector is zero.");
        matrix<T> z = y * (1.0 / norm);

        for (int iter = 0; iter < max_iter; ++iter) {
            matrix<T> B = A - matrix<T>::eye(n) * sigma;

            // Регуляризация для избежания вырожденности
            B = B + matrix<T>::eye(n) * 1e-6;

            try {
                matrix<T> y_new = solve_system(B, z);
                norm = y_new.norm();
                if (norm < 1e-10) break;
                z = y_new * (1.0 / norm);

                // Вычисление μ_i с новым критерием
                std::vector<T> mu_values;
                for (int i = 0; i < n; ++i) {
                    if (std::abs(y_new[i][0]) > delta) {
                        mu_values.push_back(z[i][0] / y_new[i][0]);
                    }
                }

                if (mu_values.empty())
                    throw std::runtime_error("All components below threshold.");

                T avg_mu = std::accumulate(mu_values.begin(), mu_values.end(), T(0)) / mu_values.size();
                T sigma_new = sigma + avg_mu;

                // Критерий сходимости
                if (std::abs(sigma_new - sigma) < epsilon && (y_new - z * norm).norm() < epsilon) {
                    return { sigma_new, z };
                }

                sigma = sigma_new;
                y = y_new;

            }
            catch (const std::exception& e) {
                throw std::runtime_error("Solver error: " + std::string(e.what()));
            }
        }
        throw std::runtime_error("Method did not converge.");
    }

    // 4. Поиск всех собственных пар с улучшенной стратегией
    template<typename T>
    std::vector<std::pair<T, matrix<T>>> find_all_eigen_pairs(const matrix<T>& A, double epsilon = 1e-6) {
        std::vector<std::pair<T, matrix<T>>> eigen_pairs;
        int n = A.getrow();
        std::vector<T> shifts = { /* Диагональные элементы как начальные сдвиги */
            A[0][0], A[1][1], A[2][2]
        };

        for (T sigma : shifts) {
            // Генерация случайного вектора с нормальным распределением
            matrix<T> vec = matrix<T>::random(n, 1, -1.0, 1.0);
            try {
                auto pair = inverse_power_method(A, sigma, vec, 1e-8, epsilon);

                // Проверка уникальности с относительным допуском
                bool unique = true;
                for (const auto& ep : eigen_pairs) {
                    if (std::abs((ep.first - pair.first) / pair.first) < 0.01) {
                        unique = false;
                        break;
                    }
                }
                if (unique) eigen_pairs.push_back(pair);
            }
            catch (...) {
                // Пропуск неудачных попыток
            }
        }
        return eigen_pairs;
    }
#endif
}

// //
//  logic of the apply Method To Elements method using SFINAE to check for the presence of a 
// method with a parameter for a matrix element and to call it if it exists

//the main template/
template <typename T, typename Param, typename = void>
struct HasMethodWithParam : std::false_type {};

//Now the compiler is thinking whether to use the main template, or somewhere there is a separate specialization for such a case.
template <typename T, typename Param>
struct HasMethodWithParam<T, Param, std::void_t<decltype(std::declval<T>().output_mode_set(std::declval<Param>()))>> : std::true_type {};
// //


template<typename T> std::ostream& operator<<(std::ostream& out, const matrix<T>& plnm);
template<typename T> std::istream& operator>>(std::istream& in, matrix<T>& plnm);
template <typename T> class matrix
{
private:
    //output_mode out_mode 
    output_mode out_mode = output_mode::FULL;

    //an array with T elements
    T* ptr;
    uint64_t colsize;
    uint64_t rowsize;

    //re-allocation of memory(does not save old values)
    void allocateMemory();

public:
     
    //CONSTRUCTORS\DESTRUCTORS

    //the copying constructor
    matrix(const matrix<T>& mtrx);
    //square matrix constructor
    matrix(uint64_t size_diag) :matrix(size_diag, size_diag) {}
    //constructor with dimensions
    matrix(uint64_t colsize, uint64_t rowsize);
    //the default constructor
    matrix();
    ~matrix();//destructor


    //DATA ACCESS
        
    uint64_t getcol()const { return colsize; }
    uint64_t getrow()const { return rowsize; }
    void setcol(uint64_t colsize) { this->colsize = colsize; this->allocateMemory(); }
    void setrow(uint64_t rowsize) { this->rowsize = rowsize; this->allocateMemory(); }
    matrix<T>& set_output_mode(output_mode mode) { this->out_mode = mode; std::cout << (int)out_mode; return *this; }
    enum class output_mode get_output_mode()const { return (this->out_mode); }
    //to index1 row access operator
    T* operator[](const uint64_t index1) const { return ptr + index1 * rowsize; }


    //ARITHMETIC OPERATORS

     // the unary operator returns a matrix with inverse (multiplied by minus 1 ) elements
    matrix<T> operator-() const;
    // binary matrix addition
    matrix<T> operator+(const matrix<T>& other) const;
    // binary matrix subtraction
    matrix<T> operator-(const matrix<T>& other) const;
    // binary matrix multiplication
    matrix<T> operator*(const matrix<T>& other) const; 
    
    // binary matrix division(if not a singular matrix on the left)
    matrix<T> operator/(const matrix<T>& other) const;
    // binary matrix multiplication by  an element from the field
    matrix<T> operator*(const T& other) const;
    // assignment operator overload
    matrix<T>& operator=(const matrix<T>& other); 


    // I/O OPERATIONS
     
    // overloading the output operator
    template<typename T>friend std::ostream& operator<<<>(std::ostream& out, const matrix<T>& p);
    // overloading the input operator
    template<typename T>friend std::istream& operator>><>(std::istream& in, matrix<T>& p);


    //SPECIAL METHODS
    
    //matrix<T> to_upper_triangular() const;

    //return of the upper triangular matrix after transformations
    matrix to_uptrng()const;
    //bringing to the upper triangular view together with the "other" matrix
    matrix to_uptrng(matrix<T>& other)const;
    //matrix transposition
    matrix transpose()const;
    // finding the determinant if there is one
    T determinant() const;
    //return of the square matrix from 1 to the lower diagonal
    //  <-----S---->
    //  0 0  ..  0 1     |
    //  0 0  ..  1 0     |
    //  . ..   ...  .. ..     S                    
    //  0 1  ..  0 0     |
    //  1 0  ..  0 0     |
    //
    matrix sqprediag(const uint64_t S)const;
    //return of the inverse matrix  
    matrix inverse_M()const;

    matrix<T> cholesky() const;
    //replacement by a matrix of zero T elements
    static matrix<T> zeros(uint64_t colsize, uint64_t rowsize);
    //replacement by a matrix of single T elements
    static matrix<T> ones(uint64_t colsize, uint64_t rowsize);
    //The identity matrix
    //  <-----S---->
    //  1 0  ..  0 0     |
    //  0 1  ..  0 0     |
    //  . ..   ...  .. ..     S                    
    //  0 0  ..  1 0     |
    //  0 0  ..  0 1     |
    static matrix<T> eye(uint64_t S);

    static matrix<T> random(size_t rows, size_t cols, T min, T max) {
        matrix<T> m(rows, cols);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dist(min, max);

        for (size_t i = 0; i < rows * cols; ++i) {
            m.ptr[i] = dist(gen);
        }
        return m;
    }

    static matrix<T> randomDiagonal(size_t n, T min, T max) {
        matrix<T> m(n, n);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dist(min, max);

        for (size_t i = 0; i < n; ++i) {
            m[i][i] = dist(gen);
        }
        return m;
    }

    // Method for calculating the maximum eigenvalue
    T max_eigenvalue(double epsilon = 1e-6, int max_iter = 1000) const;
    // Method for calculating the maximum eigenvalue with initial vector
    T max_eigenvalue(const matrix<T>& initial_vec, double epsilon = 1e-6, int max_iter = 1000) const;

    // A method for rearranging rows of a matrix
    void swap_rows(size_t i, size_t j) {
        for (size_t col = 0; col < rowsize; ++col) {
            std::swap((*this)[i][col], (*this)[j][col]);
        }
    }

    LUResult<T> LUP() const;


    //applyMethodToElements function using SFINAE to check for a method with a parameter and call it
    //for example:
    //Matrix<ExampleClass> matrix;
    //matrix.applyMethodToElements(&ExampleClass::output_mode_set, 2);
    //a method that applies a method( output_mode_set(M) ) of class T to each element of the matrix
    //(Not tested. Everything worked in the simplified code) 
    template <typename T, typename ReturnType, typename Param>
    std::enable_if_t<HasMethodWithParam<T, Param>::value> applyMethodToElements(ReturnType(T::* method)(Param), Param param)const;



    #if 1
    // Метод нормы для вектора
    
    T norm() const {
        if (rowsize != 1 && colsize != 1) {
            throw std::runtime_error("Norm is defined only for vectors.");
        }
        T sum = 0;
        for (uint64_t i = 0; i < colsize * rowsize; ++i) {
            sum += ptr[i] * ptr[i];
        }
        return std::sqrt(sum);
    }
    #endif
};
#include "../_cpp_realisation_file/matrix.cpp"
