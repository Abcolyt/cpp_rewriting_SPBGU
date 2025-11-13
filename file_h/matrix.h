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
#include <complex>

//#include "../file_h/complex.h"

extern enum class output_mode;
template<typename T> class matrix;
template<typename T>
struct is_matrix : std::false_type {};
template<typename T>
struct is_matrix<matrix<T>> : std::true_type {};

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
    template<typename T>std::pair<T, matrix<T>> power_method_max(const matrix<T>& A, const matrix<T>& Vec0, double epsilon = 1e-6, int max_iter = 1000);

    template<typename T>std::pair<T, matrix<T>> power_method_max(const matrix<T>& A, double epsilon = 1e-6, int max_iter = 1000) {
        return matrixfunction::power_method_max(A, matrix<T>::ones(A.getcol(), 1), epsilon, max_iter);
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
    
////Hessenberg form
    
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
        // Вычисление mu = 2 / (v^T * v) !
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

    template <typename T>
    matrix<T> hessenberg_upper_form(matrix<T> A) {
        int n = A.getcol();
        if (n != A.getrow()) throw std::invalid_argument("Matrix must be square");

        for (int k = 0; k < n - 2; ++k) {
            matrix<T> H = get_the_Householder_matrix_for_reduction_to_the_upper_Hessenberg_Matrix(A, k);
            A = H * A * H; // <= H==H^-1
        }
        return A;
    }

////

//// QR Decomposition
    template <typename T>
    QRResult<T> qr_decomposition(const matrix<T>& A) {
        int m = A.getrow();
        int n = A.getcol();
        QRResult<T> result;
        result.Q = matrix<T>::eye(m); 
        result.R = A;                

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
    //2*2 matrix
    template<typename T>
    std::pair<std::complex<T>, std::complex<T>> compute_2x2_eigenvalues(const matrix<T>& A) {
        T a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1];
        T trace = a + d;
        T det = a * d - b * c;
        T discriminant = trace * trace - 4 * det;

        if (discriminant >= 0) {
            T sqrt_disc = std::sqrt(discriminant);
            return { std::complex<T>((trace + sqrt_disc) / 2, 0),
                    std::complex<T>((trace - sqrt_disc) / 2, 0) };
        }
        else {
            T real = trace / 2;
            T imag = std::sqrt(-discriminant) / 2;
            return { std::complex<T>(real, imag),
                    std::complex<T>(real, -imag) };
        }
    }
    // a lot of shift

    //n*n

    template <typename T>
    std::vector<std::complex<double>> compute_eigenvalues(const matrix<T>& A, double eps = 1e-9) {
        auto H = matrixfunction::sanitize_zeros(matrixfunction::hessenberg_upper_form(A), (T)1e-10);
        std::vector<std::complex<double>> eigenvalues;

        while (H.getrow() >= 3 && H.getcol() >= 3) {
            const uint64_t N = H.getcol() - 1;

            if (std::abs(H[N][N - 1]) <= eps) {
                eigenvalues.push_back(H[N][N]);
                H = H.submatrix(0, 0, N, N);
                continue;
            }

            auto qrH = H.qr();
            H = matrixfunction::sanitize_zeros(qrH.R * qrH.Q, (T)eps);
        }

        if (H.getrow() >= 2 && H.getcol() >= 2) {
            auto last_evs = matrixfunction::compute_2x2_eigenvalues(H);
            eigenvalues.push_back(last_evs.first);
            eigenvalues.push_back(last_evs.second);
        }
        else if (H.getrow() == 1 && H.getcol() == 1) {
            eigenvalues.push_back(H[0][0]);
        }

        return eigenvalues;
    }

    template <typename T>std::vector<std::complex<double>> compute_eigenvalues_3_qr(const matrix<T>& A, double eps = 1e-12) {
#define second_convergence false
        auto H = matrixfunction::sanitize_zeros(matrixfunction::hessenberg_upper_form(A), static_cast<T>(1e-12));
        std::vector<std::complex<double>> eigenvalues;

        // Переменные для отслеживания состояния сходимости
        T prev_mu = T(0);
        bool first_iteration = true;
        uint64_t current_size = H.getrow();

        while (H.getrow() >= 3 && H.getcol() >= 3) {
            const uint64_t n = H.getrow();
            const uint64_t N = n - 1;

            if (n != current_size) {
                first_iteration = true;
                current_size = n;
            }

            if (std::abs(H[N][N - 1]) <= eps) {
                eigenvalues.push_back(H[N][N]);
                H = H.submatrix(0, 0, N, N);
                //std::cout <<"std::abs(H[N][N - 1]) <= eps"<< std::abs(H[N][N - 1])  << "\n";
                continue;
            }
#if second_convergence == true
            if (!first_iteration && std::abs(H[N][N] - prev_mu) < (1.0 / 3.0) * std::abs(prev_mu)) {
                eigenvalues.push_back(H[N][N]);
                H = H.submatrix(0, 0, N, N);
                first_iteration = true;
                //std::cout << "!first_iteration && std::abs(H[N][N] - prev_mu) < (1.0 / 3.0) * std::abs(prev_mu):\nabs:" << std::abs(H[N][N - 1]) << "\n abs(prey_mu)"<< std::abs(prev_mu)<<"\n";
                continue;
            }
#endif
            //  сдвиг 
            T mu = H[N][N];
            prev_mu = mu;  
            first_iteration = false;

            matrix<T> H_shifted = H - matrix<double>::eye(H.getrow()) * mu;
            //for (uint64_t i = 0; i < n; ++i) {
            //    H_shifted[i][i] -= mu;  // H - mu*I
            //}

            auto qrH = H_shifted.qr();      
            H = qrH.R * qrH.Q;  
            //H=H + matrix<double>::eye(H.getrow()) * mu;
            for (uint64_t i = 0; i < n; ++i) {
                H[i][i] += mu;              // + mu*I
            }

            H = matrixfunction::sanitize_zeros(H, static_cast<T>(eps));
        }

        if (H.getrow() == 2 && H.getcol() == 2) {
            auto last_evs = matrixfunction::compute_2x2_eigenvalues(H);
            eigenvalues.push_back(last_evs.first);
            eigenvalues.push_back(last_evs.second);
        }
        else if (H.getrow() == 1 && H.getcol() == 1) {
            eigenvalues.push_back(H[0][0]);
        }

        return eigenvalues;
    }


#if 0
    template <typename T>
    std::vector<T> qr_algorithm_with_shifts(matrix<T>& H, T epsilon = 1e-6, int max_iter = 1000) {
        int n = H.getcol();
        std::vector<T> eigenvalues;

        while (n > 0) {
            int iter = 0;
            T prev_shift = H[n - 1][n - 1];
            bool converged = false;

            while (iter < max_iter && !converged) {
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

                matrix<T> I = matrix<T>::eye(n);
                QRResult<T> qr = qr_decomposition(H - (I* shift));
                H = qr.R * qr.Q +  I* shift;

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

            if (!converged) {
                eigenvalues.push_back(H[n - 1][n - 1]);
                n--;
                H = H.submatrix(0, 0, n, n);
            }
        }

        return eigenvalues;
    }
#else
    template <typename T>std::vector<T> qr_algorithm_with_shifts(matrix<T>& H, T epsilon = 1e-6, int max_iter = 1000) {
        while (H.getrow() >= 2 && H.getcol())
            compute_2x2_eigenvalues(H.submatrix(H.getrow() - 2, H.getcol() - 2, H.getrow(), H.getcol()));
                return 0;
    }
#endif
    //
   

////

////scalar product of vectors
    template <typename T>
    T dot_product(const matrix<T>& a, const matrix<T>& b) {

        const bool a_is_col = (a.getcol() > 1 && a.getrow() == 1);
        const bool a_is_row = (a.getrow() > 1 && a.getcol() == 1);
        const bool b_is_col = (b.getcol() > 1 && b.getrow() == 1);
        const bool b_is_row = (b.getrow() > 1 && b.getcol() == 1);

        if (!(a_is_col || a_is_row) || !(b_is_col || b_is_row)) {
            throw std::invalid_argument("Both arguments must be vectors");
        }

        const uint64_t a_len = a_is_col ? a.getcol() : a.getrow();
        const uint64_t b_len = b_is_col ? b.getcol() : b.getrow();

        if (a_len != b_len) {
            throw std::invalid_argument("Vectors must have the same length");
        }

        T result = 0;
        for (uint64_t i = 0; i < a_len; ++i) {
            const T a_val = a_is_col ? a[0][i] : a[i][0]; // Для столбца [i][0], для строки [0][i]
            const T b_val = b_is_col ? b[0][i] : b[i][0];
            result += a_val * b_val;
        }

        return result;
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
    //2D Constructor with initialization list
    matrix(std::initializer_list<std::initializer_list<T>> init);
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
    //move operator
    matrix(matrix<T>&& other) noexcept
        : ptr(other.ptr), colsize(other.colsize), rowsize(other.rowsize), out_mode(other.out_mode) {
        other.ptr = nullptr;
        other.colsize = 0;
        other.rowsize = 0;
    }
    ~matrix();//destructor


    //DATA ACCESS
        
    uint64_t getcol()const { return colsize; }
    uint64_t getrow()const { return rowsize; }
    void setcol(uint64_t colsize) { this->colsize = colsize; this->allocateMemory(); }
    void setrow(uint64_t rowsize) { this->rowsize = rowsize; this->allocateMemory(); }
    matrix<T>& set_output_mode(output_mode mode) { this->out_mode = mode; return *this; }
    enum class output_mode get_output_mode()const { return (this->out_mode); }
    //to index row access operator
    T* operator[](const uint64_t row_index) const { return ptr + row_index * rowsize; }
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
    //r-value moving
    matrix<T>& operator=(matrix<T>&& other) noexcept {
        if (this != &other) {
            delete[] ptr;
            ptr = other.ptr;
            colsize = other.colsize;
            rowsize = other.rowsize;
            out_mode = other.out_mode;
            other.ptr = nullptr;
            other.colsize = 0;
            other.rowsize = 0;
        }
        return *this;
    }


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
    matrix<T> operator/(const T& other) const;

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
    //  <---rows--->
    //  1 1  ..  1 1     |
    //  1 1  ..  1 1     |
    //  . ..   ...  .. ..    cols                   
    //  1 1  ..  1 1     |
    //  1 1  ..  1 1     |
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
    static matrix<T> random(size_t rows, size_t cols, T min, T max);

    //Random matrix
    // min < R < max
    //  <---rows--->
    //  R 0  ..  0 0     |
    //  0 R  ..  0 0     |
    //  . ..   ...  .. ..     cols                    
    //  0 0  ..  R 0     |
    //  0 0  ..  0 R     |
    static matrix<T> randomDiagonal(size_t n, T min, T max);

    // Constructed from vector x of length n, with m columns (if m=0 then m=n)
    // Matrix dimensions: n rows x m columns
    // Element at row i, column j: x[i]^j  [increasing powers]
    //  <--- m columns --->
    //  x0^0  x0^1  x0^2  ...  x0^(m-1)   | row 0
    //  x1^0  x1^1  x1^2  ...  x1^(m-1)   | row 1
    //  ...                                | ...
    //  xn^0  xn^1  xn^2  ...  xn^(m-1)   | row n-1
    static matrix<T> vander(const std::vector<T>& x, uint64_t m=0 );
    
    // (A'A)^-1 * A'
    matrix<T> left_pseudo_reverse()const;
    // A'*( A*A')^-1
    matrix<T> right_pseudo_reverse()const;

    // General pseudo-inverse matrix (automatic selection)
    matrix<T> pseudo_inverse()const;

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
    
    //
    //mat:
    //cols: 4, rows: 4
    // 68 -94  73 -32
    // -5  -7  32  39
    //-69 -79 -91  43
    // 81  40   6  49
    //sub:
    //cols: 3, rows: 3
    // -94  73 -32
    // -7  32  39
    // -79 -91  43
    //
    matrix<T> submatrix(uint64_t start_row, uint64_t start_col, uint64_t rows, uint64_t cols) const;

    std::vector<std::complex<T>> eigenvalues_qr_double_shift(double epsilon) const;
    
    //LOGIC OPERATIONS
    //Is a vector a row or a vector a column
    bool is_vector() {
        return (this->colsize == 1) || (this->rowsize == 1);
    }
    bool is_vertical_vector() {
        return (this->rowsize == 1);
    }
    bool is_gorizontal_vector() {
        return (this->colsize == 1);
    }
    bool operator==(const matrix<T>& other) const {
        if (this->getcol() == other.getcol() && this->getrow() == other.getrow()) {

        for (size_t i = 0; i < this->getcol(); i++)
        {
            for (size_t j = 0; j < this->getrow(); j++)
            {   
                if ((*this)[j][i] != other[j][i]) {
                    return false;
                }
            }
        }
        return true;
        }
        else {
            return false;
        }
    }
};

template<typename T, typename Scalar>
auto operator*(const Scalar& scalar, const matrix<T>& mat)
-> std::enable_if_t<std::is_arithmetic_v<Scalar>, matrix<T>>
{
    return mat * static_cast<T>(scalar);
}
#include "../_cpp_realisation_file/matrix.cpp"
