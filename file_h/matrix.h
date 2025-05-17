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

//LUP system solver
   
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

    //  Решение системы 
    template<typename T>
    matrix<T> solve_system(const matrix<T>& A, const matrix<T>& b) {
        auto lup = A.LUP();
        matrix<T> Pb = lup.P * b;
        matrix<T> y = forward_substitution(lup.L, Pb);
        matrix<T> x = backward_substitution(lup.U, y);
        return x;
    }

#if 1
    template<typename T>std::vector<T> generatePoints_equally_sufficient_(int k, T a, T b) {
        T step = (b - a) / k;
        std::vector<T> points;
        for (int i = 0; i < k; ++i) {
            T x = a + i * step;
            points.emplace_back(x) ;
        }
        return points;
    }

    template <typename T>
    std::vector<std::pair<T, matrix<T>>> inverse_power_method_with_shifts(
        matrix<T> A,
        const std::vector<T>& initial_shifts,
        double epsilon = 1e-6,
        double delta = 1e-8,
        int max_iter = 10000
    ) {
        std::vector<std::pair<T, matrix<T>>> eigen_pairs;
        for (const T& sigma0 : initial_shifts) {
            ////1.1
            matrix<T> z = matrix<T>::random(A.getcol(), 1, -1.0, 1.0);
            //z = matrix<T>::ones(A.getcol(), 1); 
            z = z * (1.0 / z.norm());
            //std::cout << z << '\n';
            ////1.2
            T sigma = sigma0;
            T sigma_prev = 0;
            ////
            bool converged = false;

            
            for (int iter = 0; iter < max_iter; ++iter) {
                //std::cout << "sigma0:" << sigma0 << " iter:" << iter << "abs eps:" << std::abs(sigma - sigma_prev) << "\n";
                
                ////2.
               // std::cout << "eye:" << (matrix<T>::eye(A.getcol()) * sigma).set_output_mode(output_mode::ABBREVIATED) << "\n";
               // std::cout << "A:" << ((A*1).set_output_mode(output_mode::ABBREVIATED)) << "\nA- E*sigma" << (A - matrix<T>::eye(A.getcol()) * sigma).set_output_mode(output_mode::ABBREVIATED) << "\n";
                matrix<T> A_shifted = A - matrix<T>::eye(A.getcol()) * sigma;
                //std::cout << "det:" << std::abs(A_shifted.determinant())<<"\n";
                if (std::abs(A_shifted.determinant()) < epsilon) {
                    //std::cerr << "Warning: Matrix (A - " << sigma << "*I) is singular. Skipping this shift.\n";
                    converged = true;
                   // eigen_pairs.emplace_back(sigma, z);
                    break;
                }
                matrix<T> y = solve_system(A_shifted, z);//Y_ k+1 = y_(iter) from z_(iter-1)
                //std::cout << "sigma:" << sigma << "\n";//"\ndet:" << y.set_output_mode(output_mode::ABBREVIATED) << "\n";
                ////
                
                //T mu = (z.transpose() * y)[0][0];
                //std::cout << "(z.transpose()).getcol():" << (z.transpose()).getcol() << "(z.transpose()).getrow():" << (z.transpose()).getrow() << "\n" << (z.transpose()).set_output_mode(output_mode::FULL) << "\n";
                //std::cout << "(y).getcol():" << (y).getcol() << "(y).getrow():" << (y).getrow() << "\n" << (y).set_output_mode(output_mode::ABBREVIATED) << "\n";
                //std::cout<<"(z.transpose() * y).getcol():" << (z.transpose() * y).getcol() << "  (z.transpose() * y).getrow():" << (z.transpose() * y).getrow()<<"\n"<< (z.transpose() * y).set_output_mode(output_mode::ABBREVIATED) << "\n";
                
                //// 3.
                T mu = 0;
#if 1
                int count = 0;
                for (uint64_t i = 0; i < y.getcol(); ++i) {
                    //std::cout << "std::abs(y["<<i<<"][0]):" << std::abs(y[i][0])<<" z[i][0] / y[i][0] : "<< z[i][0] / y[i][0] << "\n";
                    //std::cout << "std::abs(y[" << i << "][0]):" << std::abs(y[i][0]) << "delta:" << delta<<">=:"<< (std::abs(y[i][0]) >= delta) << "\n";
                    if ((std::abs(y[i][0]) >= delta)) { // without 
                        mu += z[i][0] / y[i][0];
                    count++;
                    }
                }
                //count = y.getcol();
                //std::cout << "count:" << count << "\n";
                if (count == 0) break; 
                mu = mu / count; // Среднее μ
                //..std::cout << "mu:" << mu << "\n";

#else
                mu = (z.transpose() * y)[0][0];
                if (std::abs(mu) < epsilon) {
                    std::cout << "mu is too small. Breaking iteration.\n";
                    break;
                }

                mu = 1 /mu;
#endif
                sigma = sigma +  mu;
                //std::cout << "sigma:" << sigma << "\n";
                ////converged?

                if (std::abs(sigma - sigma_prev) < epsilon) {
                    converged = true;
                    break;
                }
                
                sigma_prev = sigma;

                //z_(iter) from y_(iter)
                T y_norm = y.norm();
                if (y_norm < epsilon) break;
                z = y * (1.0 / y_norm);
            }





            if (converged) {
                eigen_pairs.emplace_back(sigma, z);
            }
            else {
                std::cout << "Warning: Convergence not achieved for shift " << sigma <<" "<<z << ".\n";
            }
        }
        return eigen_pairs;
    }

#else 


    template <typename T>
    std::vector<std::pair<T, matrix<T>>> inverse_power_method_with_shifts(
        const matrix<T>& A,
        const std::vector<T>& initial_shifts,
        double epsilon = 1e-6,
        int max_iter = 1000
    ) {
        std::vector<std::pair<T, matrix<T>>> eigen_pairs;
        for (const T& sigma : initial_shifts) {
            matrix<T> z = matrix<T>::random(A.getcol(), 1, -1.0, 1.0);
            z = z * (1.0 / z.norm());

            T lambda_prev = 0;
            T lambda = 0;
            bool converged = false;

            matrix<T> A_shifted = A - matrix<T>::eye(A.getcol()) * sigma;
            T det = A_shifted.determinant();
            if (std::abs(det) < epsilon) {
                std::cerr << "Warning: Matrix (A - " << sigma << "I) is singular. Skipping this shift.\n";
                continue;
            }

            for (int iter = 0; iter < max_iter; ++iter) {
                matrix<T> y = solve_system(A_shifted, z);
                T y_norm = y.norm();
                if (y_norm < epsilon) break;

                T mu = (z.transpose() * y)[0][0];
                std::cout << "(z.transpose()).getcol():" << (z.transpose()).getcol() << "(z.transpose()).getrow():" << (z.transpose()).getrow() << "\n" << (z.transpose()).set_output_mode(output_mode::FULL) << "\n";
                std::cout << "(y).getcol():" << (y).getcol() << "(y).getrow():" << (y).getrow() << "\n" << (y).set_output_mode(output_mode::ABBREVIATED) << "\n";
                std::cout << "(z.transpose() * y).getcol():" << (z.transpose() * y).getcol() << "  (z.transpose() * y).getrow():" << (z.transpose() * y).getrow() << "\n" << (z.transpose() * y).set_output_mode(output_mode::ABBREVIATED) << "\n";

                // Вычисление μ_i = z_i^{(k-1)} / y_i^{(k)} для ненулевых компонент
                T sum_mu = 0;
                int count = 0;
                for (uint64_t i = 0; i < y.getcol() ; ++i) {
                    if (std::abs(y[i][0]) > delta) { // Игнорируем нули
                        sum_mu += z[i][0] / y[i][0];
                        count++;
                    }
                }
                if (count == 0) break; // Все компоненты близки к нулю

                T mu = sum_mu / count; // Среднее μ

                if (std::abs(lambda - lambda_prev) < epsilon) {
                    converged = true;
                    break;
                }
                lambda_prev = lambda;

                z = y * (1.0 / y_norm);//Z_ k+1
            }
            if (converged) {
                eigen_pairs.emplace_back(lambda, z);
            }
            else {
                std::cerr << "Warning: Convergence not achieved for shift " << sigma << ".\n";
            }
        }
        return eigen_pairs;
    }

#endif
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
    matrix<T>& set_output_mode(output_mode mode) { this->out_mode = mode; return *this; }
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
    //matrix<ExampleClass> matrix;
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
