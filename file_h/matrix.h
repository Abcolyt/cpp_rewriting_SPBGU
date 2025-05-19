#pragma once
#include<iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <random>
#include <limits>
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
template<typename T>
struct QRResult {
    matrix<T> Q;
    matrix<T> R;
};
namespace matrixfunction {
//// POWER METHOD
    template<typename T>T power_method(const matrix<T>& A, const matrix<T>& Vec0, double epsilon = 1e-6, int max_iter = 1000);

    template<typename T>T power_method(const matrix<T>& A, double epsilon = 1e-6, int max_iter = 1000) {
        return matrixfunction::power_method(A, matrix<T>::ones(A.getcol(), 1), epsilon, max_iter);
    }
////
    template <typename T>bool is_equal(const matrix<T>& a, const matrix<T>& b, T epsilon = static_cast<T>(1e-6));

    template <typename T>matrix<int> elementwise_equal_matrix(const matrix<T>& a, const matrix<T>& b, T epsilon = static_cast<T>(1e-6));

    template<typename T>
    matrix<T> sanitize_zeros(matrix<T> m, const T eps = std::numeric_limits<T>::epsilon() * 10) {
        

        for (uint64_t i = 0; i < m.getrow(); ++i) {
            for (uint64_t j = 0; j < m.getcol(); ++j) {
                if (std::abs(m[i][j]) < eps) {
                    m[i][j] = T(0);
                }
            }
        }
        return m;
    }

////LUP system solver
    template<typename T>
    matrix<T> forward_substitution(const matrix<T>& L, const matrix<T>& b);
    template<typename T>
    matrix<T> backward_substitution(const matrix<T>& U, const matrix<T>& y);
    template<typename T>
    matrix<T> solve_system(const matrix<T>& A, const matrix<T>& b);
////
    
//// QR Decomposition
    template <typename T>
    QRResult<T> qr_decomposition(const matrix<T>& A) {
        int m = A.getrow();
        int n = A.getcol();
        QRResult<T> result;
        result.Q = matrix<T>::eye(m); // Инициализация Q как единичной матрицы
        result.R = A;                 // Копирование исходной матрицы в R

        for (int j = 0; j < n; ++j) {
            for (int i = j + 1; i < m; ++i) {
                // Вычисление вращения Гивенса для обнуления R[i][j]
                T a = result.R[j][j];
                T b = result.R[i][j];
                if (std::abs(b) < 1e-10) continue;

                T r = std::hypot(a, b);
                T c = a / r;
                T s = -b / r;

                // Применение вращения к R
                for (int k = j; k < n; ++k) {
                    T temp = c * result.R[j][k] - s * result.R[i][k];
                    result.R[i][k] = s * result.R[j][k] + c * result.R[i][k];
                    result.R[j][k] = temp;
                }

                // Применение вращения к Q
                for (int k = 0; k < m; ++k) {
                    T temp = c * result.Q[k][j] - s * result.Q[k][i];
                    result.Q[k][i] = s * result.Q[k][j] + c * result.Q[k][i];
                    result.Q[k][j] = temp;
                }
            }
        }

        return result;
    }
////

////Eigenvalues and vectors
    //Simplification column vector or row vector
    template <typename T>matrix<T> simplify_eigenvector(const matrix<T>& vec, T epsilon = 1e-6);
    
    //
    template<typename T>std::vector<T> generatePoints_equally_sufficient_(int k, T a, T b) {
        T step = (b - a) / k;
        std::vector<T> points;
        for (int i = 0; i < k; ++i) {
            T x = a + i * step;
            points.emplace_back(x) ;
        }
        return points;
    }
    //a lot of shifts
    template <typename T>
    std::vector<std::pair<T, matrix<T>>> inverse_power_method_with_shifts(
        matrix<T> A,const std::vector<T>& initial_shifts,
        double epsilon = 1e-6,double delta = 1e-8,int max_iter = 10000
    );

    //one shift
    template <typename T>
    std::pair<T, matrix<T>> inverse_power_method_with_shift(
        matrix<T> A,T sigma0,const matrix<T>& vec0,
        double epsilon = 1e-6,double delta = 1e-8,int max_iter = 10000
    );

    // a lot of shift
    template <typename T>
    std::vector<T> qr_algorithm_with_shifts(matrix<T>& H, T epsilon = 1e-6, int max_iter = 1000) {
        int n = H.getcol();
        std::vector<T> eigenvalues;

        while (n > 0) {
            int iter = 0;
            T prev_shift = H[n - 1][n - 1];
            bool converged = false;

            while (iter < max_iter && !converged) {
                // Выбор сдвига (Wilkinson для блоков 2x2)
                T shift;
                bool is_2x2_block = (n >= 2 && std::abs(H[n - 2][n - 1]) > epsilon);

                if (is_2x2_block) {
                    T a = H[n - 2][n - 2];
                    T b = H[n - 2][n - 1];
                    T c = H[n - 1][n - 2];
                    T d = H[n - 1][n - 1];
                    T delta = (a - d) / 2;
                    T sqrt_term = std::sqrt(delta * delta + b * c);
                    shift = d + delta - std::copysign(sqrt_term, delta);
                }
                else {
                    shift = H[n - 1][n - 1];
                }

                // QR-разложение со сдвигом
                matrix<T> I = matrix<T>::eye(n);
                QRResult<T> qr = qr_decomposition(H - (I* shift));
                H = qr.R * qr.Q +  I* shift;

                // Проверка сходимости
                T current_shift = H[n - 1][n - 1];
                T subdiag = (n > 1) ? std::abs(H[n - 1][n - 2]) : 0;

                if (n == 1) {
                    eigenvalues.push_back(H[0][0]);
                    converged = true;
                }
                else if (subdiag < epsilon && std::abs(current_shift - prev_shift) < 0.33 * std::abs(prev_shift)) {
                    eigenvalues.push_back(current_shift);
                    n--;
                    H = H.submatrix(0, 0, n, n);
                    converged = true;
                }

                prev_shift = current_shift;
                iter++;
            }

            // Если не сошлось за max_iter, принудительно уменьшаем размерность
            if (!converged) {
                eigenvalues.push_back(H[n - 1][n - 1]);
                n--;
                H = H.submatrix(0, 0, n, n);
            }
        }

        return eigenvalues;
    }
////

////
    //scalar product of vectors
    template <typename T>
    T dot_product(const matrix<T>& a, const matrix<T>& b) {
        // Проверка что оба аргумента - векторы (столбцы или строки)
        const bool a_is_col = (a.getcol() > 1 && a.getrow() == 1);
        const bool a_is_row = (a.getrow() > 1 && a.getcol() == 1);
        const bool b_is_col = (b.getcol() > 1 && b.getrow() == 1);
        const bool b_is_row = (b.getrow() > 1 && b.getcol() == 1);

        if (!(a_is_col || a_is_row) || !(b_is_col || b_is_row)) {
            throw std::invalid_argument("Both arguments must be vectors");
        }

        // Определение длины векторов
        const uint64_t a_len = a_is_col ? a.getcol() : a.getrow();
        const uint64_t b_len = b_is_col ? b.getcol() : b.getrow();

        if (a_len != b_len) {
            throw std::invalid_argument("Vectors must have the same length");
        }

        // Вычисление скалярного произведения
        T result = 0;
        for (uint64_t i = 0; i < a_len; ++i) {
            const T a_val = a_is_col ? a[0][i] : a[i][0]; // Для столбца [i][0], для строки [0][i]
            const T b_val = b_is_col ? b[0][i] : b[i][0];
            result += a_val * b_val;
        }

        return result;
    }
////

////

#if 0

#else
    
    template <typename T>
    matrix<T> get_the_Householder_matrix_for_reduction_to_the_upper_Hessenberg_Matrix(const matrix<T>& A, int k) {
        int n = A.getcol();
        if (n != A.getrow()) throw std::invalid_argument("Matrix must be square");
        if (k < 0 || k >= n - 1) throw std::invalid_argument("Invalid column index");

        // Выбор подвектора x из столбца k, начиная с элемента k+1
        int m = n - k - 1;
        matrix<T> x(m, 1); // Столбец-вектор
        for (int i = 0; i < m; ++i) {
            x[i][0] = A[k + 1 + i][k];
        }
        //std::cout << x << "\n";

        // Вычисление нормы x и s
        T norm_x = x.norm();
        if (norm_x == 0) return matrix<T>::eye(n);

        T s = std::copysign(norm_x, x[0][0]);

        // Построение вектора v = x - s * e
        matrix<T> v = x;
        v[0][0] = v[0][0] - s;

        //std::cout << "v:" << v << "\n";
        // Вычисление mu = 2 / (v^T * v)
        matrix<T> vT = v.transpose();
        matrix<T> vtv = vT * v;
        if (vtv[0][0] == 0) return matrix<T>::eye(n);
        T mu = T(2) / vtv[0][0];

        //// Построение матрицы Хаусхолдера 1. H = I - mu * v * v^T
        matrix<T> H = matrix<T>::eye(n);
        matrix<T> outer = v * vT;
        outer = outer * mu;
        //std::cout <<"outer:" << outer << "\n";

        // 2.
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                H[k + 1 + i][k + 1 + j] -= outer[i][j];
            }
        }
        ////
        return H;
    }
#endif
    template <typename T>
    matrix<T> hessenberg_upper_form(matrix<T> A) {
        int n = A.getcol();
        if (n != A.getrow()) throw std::invalid_argument("Matrix must be square");

        for (int k = 0; k < n - 2; ++k) {
            matrix<T> H = get_the_Householder_matrix_for_reduction_to_the_upper_Hessenberg_Matrix(A, k);
            A = H * A * H; // Применение преобразования
        }
        return A;
    }

////
}



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
    output_mode out_mode = output_mode::ABBREVIATED;

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

    template <typename U>matrix(const matrix<U>& other) {
        colsize = other.getcol();
        rowsize = other.getrow();
        allocateMemory();
        for (uint64_t i = 0; i < colsize; ++i) {
            for (uint64_t j = 0; j < rowsize; ++j) {
                (*this)[i][j] = static_cast<T>(other[i][j]);
            }
        }
    }
    ~matrix();//destructor


    //DATA ACCESS
        
    uint64_t getcol()const { return colsize; }
    uint64_t getrow()const { return rowsize; }
    void setcol(uint64_t colsize) { this->colsize = colsize; this->allocateMemory(); }
    void setrow(uint64_t rowsize) { this->rowsize = rowsize; this->allocateMemory(); }
    matrix<T>& set_output_mode(output_mode mode) { this->out_mode = mode; return *this; }
    enum class output_mode get_output_mode()const { return (this->out_mode); }
    //to index1 row access operator
    T* operator[](const uint64_t index1) const { return ptr + index1 * rowsize; }
    matrix<T> get_column(uint64_t col) const;
    matrix<T> get_row(uint64_t row) const;

    //ARITHMETIC OPERATORS

     // the unary operator returns a matrix with inverse (multiplied by minus 1 ) elements
    matrix<T> operator-() const;
    // binary matrix addition
    matrix<T> operator+(const matrix<T>& other) const;
    // binary matrix subtraction
    matrix<T> operator-(const matrix<T>& other) const;
    // binary matrix multiplication
    matrix<T> operator*(const matrix<T>& other) const; 
    // binary matrix multiplication by  an element from the field
    matrix<T> operator*(const T& other) const;
    // assignment operator overload
    matrix<T>& operator=(const matrix<T>& other); 

    template <typename U>matrix<T>& operator=(const matrix<U>& other) {
        if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
            delete[] ptr;
            colsize = other.getcol();
            rowsize = other.getrow();
            allocateMemory();
            for (uint64_t i = 0; i < colsize; ++i) {
                for (uint64_t j = 0; j < rowsize; ++j) {
                    (*this)[i][j] = static_cast<T>(other[i][j]);
                }
            }
        }
        return *this;
    }
    // binary matrix division(if not a singular matrix on the left)
    matrix<T> operator/(const matrix<T>& other) const;


    // I/O OPERATIONS
     
    // overloading the output operator
    template<typename T>friend std::ostream& operator<<<>(std::ostream& out, const matrix<T>& p);
    // overloading the input operator
    template<typename T>friend std::istream& operator>><>(std::istream& in, matrix<T>& p);


    //SPECIAL METHODS
    // The p norm method for a matrix
    T norm(int p = 2) const;
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
    matrix sqprediag(const uint64_t S)const;
    //return of the inverse matrix  
    matrix inverse_M()const;
    //The Cholesky decomposition
    matrix<T> cholesky() const;
    //LUP decomposition
    LUResult<T> LUP() const;
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
    //Random matrix
    // min < R < max
    //  <---rows--->
    //  R R  ..  R R     |
    //  R R  ..  R R     |
    //  . ..   ...  .. ..    cols                   
    //  R R  ..  R R     |
    //  R R  ..  R R     |
    static matrix<T> random(size_t rows, size_t cols, T min, T max) {
        matrix<T> m(rows, cols);
        std::random_device rd;
        std::mt19937 gen(rd());

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            for (size_t i = 0; i < rows * cols; ++i) {
                m.ptr[i] = dist(gen);
            }
        }
        else {
            static_assert(
                std::is_same_v<T, float> ||
                std::is_same_v<T, double> ||
                std::is_same_v<T, long double>,
                "T must be float, double, or long double for real distribution"
                );
            std::uniform_real_distribution<T> dist(min, max);
            for (size_t i = 0; i < rows * cols; ++i) {
                m.ptr[i] = dist(gen);
            }
        }
        return m;
    }

    //Random matrix
    // min < R < max
    //  <---rows--->
    //  R 0  ..  0 0     |
    //  0 R  ..  0 0     |
    //  . ..   ...  .. ..     cols                    
    //  0 0  ..  R 0     |
    //  0 0  ..  0 R     |
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



    //applyMethodToElements function using SFINAE to check for a method with a parameter and call it
    //for example:
    //matrix<ExampleClass> matrix;
    //matrix.applyMethodToElements(&ExampleClass::output_mode_set, 2);
    //a method that applies a method( output_mode_set(M) ) of class T to each element of the matrix
    //(Not tested. Everything worked in the simplified code) 
    template <typename T, typename ReturnType, typename Param>
    std::enable_if_t<HasMethodWithParam<T, Param>::value> applyMethodToElements(ReturnType(T::* method)(Param), Param param)const;

    
    QRResult<T> qr() const;

    matrix<T> submatrix(uint64_t start_row, uint64_t start_col, uint64_t rows, uint64_t cols) const;

};
#include "../_cpp_realisation_file/matrix.cpp"
