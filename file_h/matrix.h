#pragma once
#include<iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <random>

extern enum class output_mode;
template<typename T> class matrix;

namespace matrixfunction {
    template<typename T>T power_method(const matrix<T>& A, const matrix<T>& Vec0, double epsilon = 1e-6, int max_iter = 1000);

    template<typename T>T power_method(const matrix<T>& A, double epsilon = 1e-6, int max_iter = 1000) {
        return matrixfunction::power_method(A, matrix<T>::ones(A.getcol(), 1), epsilon, max_iter);
    }
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
    void set_output_mode(output_mode mode) { out_mode = mode; }
    enum class output_mode get_output_mode()const { return this->out_mode; }
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
    
    matrix<T> to_upper_triangular() const;

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

    //applyMethodToElements function using SFINAE to check for a method with a parameter and call it
    //for example:
    //Matrix<ExampleClass> matrix;
    //matrix.applyMethodToElements(&ExampleClass::output_mode_set, 2);
    //a method that applies a method( output_mode_set(M) ) of class T to each element of the matrix
    //(Not tested. Everything worked in the simplified code) 
    template <typename T, typename ReturnType, typename Param>
    std::enable_if_t<HasMethodWithParam<T, Param>::value> applyMethodToElements(ReturnType(T::* method)(Param), Param param)const;
    

    matrix<T> cholesky() const;
    //replacement by a matrix of zero T elements
    static matrix<T> zeros(uint64_t colsize, uint64_t rowsize) {
        matrix<T> mat;
        mat.colsize = colsize;
        mat.rowsize = rowsize;
        mat.allocateMemory();
        return mat;
    }

    //replacement by a matrix of single T elements
    static matrix<T> ones(uint64_t colsize, uint64_t rowsize) {
        matrix<T> mat;
        mat.colsize = colsize;
        mat.rowsize = rowsize;
        mat.ptr = new T[colsize * rowsize];
        for (uint64_t i = 0; i < colsize * rowsize; i++) {
            mat.ptr[i] = 1;
        }
        return mat;
    }
   

    // Method for calculating the maximum eigenvalue
    T max_eigenvalue(double epsilon = 1e-6, int max_iter = 1000) const {
        if (colsize != rowsize) {
            throw std::invalid_argument("Matrix must be square.");
        }
        return matrixfunction::power_method(*this, matrix<T>::ones(colsize, 1), epsilon, max_iter);
    }
    // Method for calculating the maximum eigenvalue with initial vector
    T max_eigenvalue(const matrix<T>& initial_vec, double epsilon = 1e-6, int max_iter = 1000) const {
        return matrixfunction::power_method(*this, initial_vec, epsilon, max_iter);
    }

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

    static matrix<T> eye(uint64_t n) {
        matrix<T> mat(n, n);
        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < n; ++j) {
                mat[i][j] = 0; // Явное обнуление
            }
            mat[i][i] = 1; // Заполнение диагонали
        }
        return mat;
    }
    // Метод для перестановки строк матрицы
    void swap_rows(size_t i, size_t j) {
        for (size_t col = 0; col < rowsize; ++col) {
            std::swap((*this)[i][col], (*this)[j][col]);
        }
    }
    // Структура для возврата LU-разложения с матрицей перестановок
    struct LUResult {
        matrix<T> L;
        matrix<T> U;
        matrix<T> P; // Матрица перестановок
    };

    LUResult lu() const {
        LUResult result;
        const size_t n = colsize;
        result.L = matrix<T>::zeros(n, n);
        result.U = *this;
        result.P = matrix<T>::eye(n); // Единичная матрица

        for (size_t k = 0; k < n; ++k) {
            // Частичный выбор ведущего элемента
            size_t max_row = k;
            T max_val = std::abs(result.U[k][k]);
            for (size_t i = k + 1; i < n; ++i) {
                if (std::abs(result.U[i][k]) > max_val) {
                    max_val = std::abs(result.U[i][k]);
                    max_row = i;
                }
            }

            // Перестановка строк в U и P
            if (max_row != k) {
                result.U.swap_rows(k, max_row);
                result.P.swap_rows(k, max_row); // Синхронизируем P
            }

            // Заполнение L и U
            result.L[k][k] = 1;
            for (size_t i = k + 1; i < n; ++i) {
                result.L[i][k] = result.U[i][k] / result.U[k][k];
                for (size_t j = k; j < n; ++j) {
                    result.U[i][j] -= result.L[i][k] * result.U[k][j];
                }
            }
        }

        return result;
    }
};
#include "../_cpp_realisation_file/matrix.cpp"
