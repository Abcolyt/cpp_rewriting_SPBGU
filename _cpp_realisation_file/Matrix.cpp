﻿#pragma once
#include "../file_h/matrix.h"
template<typename T>matrix<T>::matrix(const matrix<T>& mtrx) 
{
    (*this).ptr = nullptr;
    (*this) = mtrx;
}
template<typename T>matrix<T>::matrix() 
{
    this->colsize = 0;
    this->rowsize = 0;
    ptr = nullptr;
}

template<typename T>matrix<T>::matrix(uint64_t colsize, uint64_t rowsize)
{
    this->colsize = colsize;
    this->rowsize = rowsize;
    ptr = new T[colsize * rowsize];
    for (uint64_t i = 0; i < colsize * rowsize; i++)
    {
        ptr[i] = 0;
    }
}
template<typename T>matrix<T>::~matrix()
{
    delete[] ptr;
}

template<typename T>matrix<T> matrix<T>::sqprediag(const uint64_t S)const
{
    matrix<T> temp(S, S);
    for (int i = S - 1; i >= 0; i--)
    {
        temp[i][S - 1 - i] = 1;
    }
    return temp;
}

template<typename T>matrix<T> matrix<T>::transpose() const {

    matrix<T> temp(rowsize, colsize);

    for (size_t i = 0; i < colsize; ++i) {
        for (size_t j = 0; j < rowsize; ++j) {
            temp[j][i] = (*this)[i][j];
        }
    }
    return temp;
}

template<typename T>matrix<T> matrix<T>::inverse_M()const
{
    if (colsize != rowsize)
    {
        throw std::invalid_argument("the number of rows in the matrix is not equal to the number of columns");
    }
    uint64_t sizeM = (*this).colsize;
    matrix<T> a = (*this), reversM, E_lower_diagonal, t1;
    E_lower_diagonal = a.sqprediag(sizeM);
    reversM = E_lower_diagonal * E_lower_diagonal;

    t1 = a.to_uptrng(reversM);
    

    t1 = ((t1 * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();
    reversM = ((reversM * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();

    t1 = t1.to_uptrng(reversM);

    t1 = ((t1 * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();
    reversM = ((reversM * E_lower_diagonal).transpose() * E_lower_diagonal).transpose();

    for (uint64_t i = 0; i < sizeM; i++)
    {
        T koef = t1[i][i];
        if (koef == 0)
        {
            throw std::invalid_argument("the matrix is not reversible");
        }
        for (uint64_t k = 0; k < sizeM; k++)
        {
            t1[i][k] = t1[i][k] / koef;
            reversM[i][k] = reversM[i][k] / koef;
        }

    }
    return reversM;
}

template<typename T>matrix<T> matrix<T>::to_uptrng(matrix<T>& other)const
{
    if ((this->colsize == this->rowsize) && (this->colsize == (other.colsize)) && ((other.colsize) == other.rowsize))
    {
        //std::cout << "std::invalid_argument(the matrix is irreducible to the triangular form)\n";



        matrix<T> temp((*this));
        for (uint64_t i = 0; i < rowsize - 1; i++) {
            for (uint64_t j = i + 1; j < rowsize; j++) {

                while (temp[i][i] == 0 && temp[j][i] == 0) {
                    j++;
                    if (j >= rowsize) {
                        throw std::invalid_argument("the matrix is irreducible to the triangular form");
                        

                    }
                }

                if (temp[i][i] == 0) {
                    for (uint64_t k = 0; k < rowsize; k++) {
                        std::swap(temp[i][k], temp[j][k]);
                        std::swap(other[i][k], other[j][k]);

                    }

                }

                T koef1 = temp[j][i] / temp[i][i];
                for (uint64_t k = 0; k < rowsize; k++) {
                    temp[j][k] = temp[j][k] - temp[i][k] * koef1;
                    other[j][k] = other[j][k] - other[i][k] * koef1;
                }

            }

        }
        return temp;
    }
    throw std::invalid_argument("[  ((this->colsize == this->rowsize) && (this->colsize == (other.colsize)) && ((other.colsize) == other.rowsize)) ]==0");
  
}

template<typename T>matrix<T> matrix<T>::to_uptrng()const
{
    matrix<T> temp((*this));

    T det;
    det = 1;
    for (uint64_t i = 0; i < rowsize - 1; i++) {
        for (uint64_t j = i + 1; j < rowsize; j++) {
            while (temp[i][i] == 0 && temp[j][i] == 0) {
                j++;
                if (j >= rowsize) {
                    std::cout << "std::invalid_argument(the matrix is irreducible to the triangular form)\n";

                    throw std::invalid_argument("the matrix is irreducible to the triangular form");


                }
            }
            if (temp[i][i] == 0) {
                for (uint64_t k = 0; k < rowsize; k++) {
                    det = det * (-1);
                }
            }

            T koef = temp[j][i] / temp[i][i];
            for (uint64_t k = i; k < rowsize; k++) {
                temp[j][k] = temp[j][k] - temp[i][k] * koef;

            }

        }

    }
    for (uint64_t i = 0; i < rowsize - 1; i++) {        
        temp[0][i] = (temp[0][i]) * det;

    }
    return temp;
}


template <typename T>T matrix<T>::determinant() const {
    if (rowsize != colsize) {
        
        throw std::invalid_argument("rowsize != colsize");
    }

    matrix<T> temp((*this).to_uptrng());

    T det;
    det = (T)(1);
    //std::cout <<"temp:\n" << temp << "\n";
    for (uint64_t i = 0; i < rowsize; i++) {
        
        det = (det *temp[i][i]);
        //std::cout << (temp[i][i]) << '\n';
    }
    //std::cout << det <<'\n';
    return det;
}

template<typename T>std::ostream& operator<<(std::ostream& out, const matrix<T>& mtrx) {
    const auto mode = mtrx.get_output_mode();
    const uint64_t cols = mtrx.getcol();
    const uint64_t rows = mtrx.getrow();

    // Определяем максимальную ширину для каждого столбца
    std::vector<size_t> col_widths(rows, 1);
    if (mode != output_mode::SHORT) {
        for (uint64_t j = 0; j < rows; ++j) { // Идем по столбцам матрицы
            size_t max_width = 0;
            for (uint64_t i = 0; i < cols; ++i) { // Идем по строкам матрицы
                std::ostringstream oss;
                oss << mtrx[i][j]; // Элемент [i][j] - строка i, столбец j
                max_width = std::max(max_width, oss.str().size());
            }
            col_widths[j] = max_width + 1; // Добавляем пробел для разделения
        }
    }

    // Вывод размеров матрицы
    if (mode != output_mode::SHORT) {
        out << "cols: " << cols << ", rows: " << rows << "\n";
    }

    // Вывод содержимого
    switch (mode) {
    case output_mode::FULL: {
        for (uint64_t i = 0; i < cols; ++i) {
            for (uint64_t j = 0; j < rows; ++j) {
                out << "[" << std::setw(2) << i << "]["
                    << std::setw(2) << j << "] = "
                    << std::setw(col_widths[j]) << mtrx[i][j] << " | ";
            }
            out << "\n";
        }
        break;
    }

    case output_mode::ABBREVIATED: {
        for (uint64_t i = 0; i < cols; ++i) {
            for (uint64_t j = 0; j < rows; ++j) {
                out << std::setw(col_widths[j]) << mtrx[i][j];
            }
            out << "\n";
        }
        break;
    }

    case output_mode::SHORT: {
        out << "Matrix[" << cols << "x" << rows << "]";
        break;
    }
    }

    return out;
}


template<typename T>std::istream& operator>>(std::istream& in, matrix<T>& mtrx)
{
    uint64_t size;
    std::cout << "Enter number of rows: ";
    in >> size;
    mtrx.setrow(size);

    std::cout << "Enter number of columns: ";
    in >> size;
    mtrx.setcol(size);
    mtrx.allocateMemory();
    for (uint64_t i = 0; i < mtrx.getrow(); i++) {
        for (uint64_t j = 0; j < mtrx.getcol(); j++) {
            std::cout << "["<<i<<"]["<<j<<"]=";
            in >> mtrx[i][j];
        }
    }
    
    return in;
}

template<typename T>matrix<T>& matrix<T>::operator=(const matrix<T>& other)
{
    if (this != &other)
    {

        delete[] ptr;
        ptr = new T[other.getcol() * other.getrow()];
        this->colsize = other.getcol();
        this->rowsize = other.getrow();

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                (*this)[i][j] = other[i][j];
            }
        }
    }

    return *this;
}
template<typename T>matrix<T> matrix<T>::operator-() const
{
    matrix<T> result(this->getcol(), this->getrow());
    

        for (uint64_t i = 0; i < this->getrow(); i++)
        {
            for (uint64_t j = 0; j < this->getcol(); j++)
            {
                result[i][j] = (-1)*((*this)[i][j]);
            }

        }

    
    return result;

}
template<typename T>matrix<T> matrix<T>::operator+(const matrix<T>& other) const
{
    matrix<T> result(other.getcol(), other.getrow());
    if ((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) {

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                result[i][j] = (other[i][j] + (*this)[i][j]);
            }

        }

    }
    else {
        throw std::invalid_argument("((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) == false");
    }
    return result;

}

template<typename T>matrix<T> matrix<T>::operator-(const matrix<T>& other) const
{
    matrix<T> result(other.getcol(), other.getrow());
    if ((other.getcol()) == ((*this).getcol()) && (other.getrow() == (*this).getrow())) {

        for (uint64_t i = 0; i < other.getrow(); i++)
        {
            for (uint64_t j = 0; j < other.getcol(); j++)
            {
                result[i][j] = ((*this)[i][j]- other[i][j]);
            }

        }

    }
    return result;

}

template<typename T>matrix<T> matrix<T>::operator*(const T& other) const
{
    matrix<T> result(this->getcol(), this->getrow());


    for (uint64_t i = 0; i < this->getrow(); i++)
    {
        for (uint64_t j = 0; j < this->getcol(); j++)
        {
            result[i][j] = other * ((*this)[i][j]);
        }

    }


    return result;
}

template<typename T>matrix<T> matrix<T>::operator*(const matrix<T>& other) const
{
    if (rowsize != other.colsize)
    {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");

    }

    matrix<T> result(colsize, other.rowsize);

    for (uint64_t i = 0; i < colsize; ++i) {
        for (uint64_t j = 0; j < other.rowsize; ++j) {
            for (uint64_t k = 0; k < rowsize; ++k) {
                result[i][j] = result[i][j] + ((*this)[i][k] * other[k][j]);
            }
        }
    }

    return result;
}

template<typename T>matrix<T> matrix<T>::operator/(const matrix<T>& other) const
{
    if (other.determinant() != 0) {
        return (*this) * (other.inverse_M());
    }
    throw std::invalid_argument("the matrix on the right is not reversible");
}

template<typename T>void matrix<T>::allocateMemory()
{
    delete[] ptr;
    ptr = new T[colsize * rowsize];
    for (uint64_t i = 0; i < colsize * rowsize; i++)
    {
        ptr[i] = 0;
    }
}

template<typename T>matrix<T> matrix<T>::cholesky() const
{
    if (colsize != rowsize) {
        throw std::invalid_argument("Matrix must be square for Cholesky decomposition.");
    }

    matrix<T> L(colsize, rowsize); // Создаем матрицу L размером n x n

    for (uint64_t i = 0; i < colsize; ++i) {
        for (uint64_t j = 0; j <= i; ++j) {
            T sum = 0;

            for (uint64_t k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                // Диагональные элементы

               //std::cout << "i=" << i << "j=" << j << " " << ((*this)[i][i] - sum) << "   ";


                L[i][j] = std::sqrt((*this)[i][i] - sum);


                //std::cout << "i=" << i << "j=" << j << " " << L[i][j] <<"   ";
            }
            else {
                // Недиагональные элементы
                L[i][j] = ((*this)[i][j] - sum) / (L[j][j]);

            }
        }
    }

    return L;
}
template<typename T>LUResult<T> matrix<T>::LUP() const {
    if (colsize != rowsize) {
        throw std::invalid_argument("LU decomposition requires square matrix");
    }

    const size_t n = colsize;
    LUResult<T> result;
    result.L = matrix<T>::zeros(n, n);
    result.U = *this; // Копируем исходную матрицу
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
            result.P.swap_rows(k, max_row);

            // Перестановка строк в L (только уже вычисленные элементы)
            if (k > 0) {
                for (size_t j = 0; j < k; ++j) {
                    std::swap(result.L[k][j], result.L[max_row][j]);
                }
            }
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



template<typename T>
 matrix<T> matrix<T>::zeros(uint64_t colsize, uint64_t rowsize) {
    matrix<T> mat;
    mat.colsize = colsize;
    mat.rowsize = rowsize;
    mat.allocateMemory();
    return mat;
}

template<typename T>
 matrix<T> matrix<T>::ones(uint64_t colsize, uint64_t rowsize) {
    matrix<T> mat;
    mat.colsize = colsize;
    mat.rowsize = rowsize;
    mat.ptr = new T[colsize * rowsize];
    for (uint64_t i = 0; i < colsize * rowsize; i++) {
        mat.ptr[i] = 1;
    }
    return mat;
}
 
 template<typename T>
 matrix<T> matrix<T>::eye(uint64_t S) {
     matrix<T> mat(S, S);
     for (uint64_t i = 0; i < S; ++i) {
         mat[i][i] = 1;
     }
     return mat;
 }


 // Method for calculating the maximum eigenvalue
 template<typename T>
 T matrix<T>::max_eigenvalue(double epsilon , int max_iter) const {
     if (colsize != rowsize) {
         throw std::invalid_argument("Matrix must be square.");
     }
     return matrixfunction::power_method(*this, matrix<T>::ones(colsize, 1), epsilon, max_iter);
 }
 // Method for calculating the maximum eigenvalue with initial vector
 template<typename T>
 T matrix<T>::max_eigenvalue(const matrix<T>& initial_vec, double epsilon , int max_iter ) const {
     return matrixfunction::power_method(*this, initial_vec, epsilon, max_iter);
 }
template<typename T>
 T matrix<T>::norm(int p) const {
     if (p < 1 && p != INT_MAX) {
         throw std::invalid_argument("Norm order must be >= 1 or 'inf'.");
     }

     T result = 0;
     const uint64_t size = rowsize * colsize;

     // L1 
     if (p == 1) {
         for (uint64_t i = 0; i < size; ++i) {
             result += std::abs(ptr[i]);
         }
     }
     // L2 
     else if (p == 2) {
         for (uint64_t i = 0; i < size; ++i) {
             result += ptr[i] * ptr[i];
         }
         result = std::sqrt(result);
     }
     // Inf norm
     else if (p == INT_MAX) {
         for (uint64_t i = 0; i < size; ++i) {
             T abs_val = std::abs(ptr[i]);
             if (abs_val > result) {
                 result = abs_val;
             }
         }
     }
     //Lp
     else {
         for (uint64_t i = 0; i < size; ++i) {
             result += std::pow(std::abs(ptr[i]), p);
         }
         result = std::pow(result, 1.0 / p);
     }

     return result;
 }

 template <typename T>
 matrix<T> matrix<T>::get_column(uint64_t col) const {
     if (col >= colsize) {  
         throw std::out_of_range("Column index out of bounds");
     }

     matrix<T> column_vec(rowsize, 1); // vec: rowsize x 1
     for (uint64_t i = 0; i < rowsize; ++i) {
         column_vec[i][0] = (*this)[i][col];
     }
     return column_vec;
 }

 template <typename T>
 matrix<T> matrix<T>::get_row(uint64_t row) const {
     if (row >= rowsize) {  // <-- Исправлено colsize → rowsize
         throw std::out_of_range("Row index out of bounds");
     }

     matrix<T> row_vec(1, colsize); // Вектор-строка 1 x colsize
     for (uint64_t j = 0; j < colsize; ++j) {
         row_vec[0][j] = (*this)[row][j];
     }
     return row_vec;
 }

 template <typename T>
 QRResult<T> matrix<T>::qr() const {
     uint64_t m = getcol();
     uint64_t n = getrow();

     matrix<T> Q = matrix<T>::zeros(m, n);
     matrix<T> R = matrix<T>::zeros(n, n);

     for (uint64_t j = 0; j < n; ++j) {
         matrix<T> v = this->get_column(j); // Используем get_column()

         for (uint64_t i = 0; i < j; ++i) {
             matrix<T> Q_i = Q.get_column(i); // Получаем i-й столбец Q
             R[i][j] = matrixfunction::dot_product(Q_i, v);   // Используем функцию dot_product

             // Вычитаем проекцию
             for (uint64_t k = 0; k < m; ++k) {
                 v[k][0] -= R[i][j] * Q_i[k][0];
             }
         }

         R[j][j] = v.norm();

         if (R[j][j] != 0) {
             for (uint64_t k = 0; k < m; ++k) {
                 Q[k][j] = v[k][0] / R[j][j];
             }
         }
     }

     return QRResult<T>{Q, R};
 }

 template <typename T>
 matrix<T> matrix<T>::submatrix(uint64_t start_row, uint64_t start_col, uint64_t rows, uint64_t cols) const {
     // Проверка корректности входных параметров
     if (start_row + rows > rowsize || start_col + cols > colsize) {
         throw std::out_of_range("Submatrix dimensions exceed matrix bounds");
     }

     matrix<T> sub(rows, cols); // Создаём новую матрицу нужного размера

     for (uint64_t i = 0; i < rows; ++i) {
         for (uint64_t j = 0; j < cols; ++j) {
             // Копируем элементы из исходной матрицы
             sub[i][j] = (*this)[start_row + i][start_col + j];
         }
     }

     return sub;
 }
namespace matrixfunction {
    template<typename T>T power_method(const matrix<T>& A, const matrix<T>& Vec0, double epsilon, int max_iter ) {

        if (A.getcol() != A.getrow()) {
            throw std::invalid_argument("Matrix must be square for power method.");
        }
        int n = A.getcol();

        T eigenvalue_prev = 0;
        T eigenvalue = 0;
        matrix<T> eigenvector = Vec0;

        for (int iter = 0; iter < max_iter; ++iter) {
            matrix<T> Ab = A * eigenvector;
#if 0
            std::cout << iter << "\n";
#endif
            // Finding the maximum modulo element
            eigenvalue = Ab[0][0];
            for (int i = 0; i < n; ++i) {
                if (std::abs(Ab[i][0]) > std::abs(eigenvalue)) {
                    eigenvalue = Ab[i][0];
                }
            }

            if (((eigenvalue) < 0 ? (eigenvalue) * (-1) : (eigenvalue)) < epsilon) {
                throw std::runtime_error("Matrix may be singular (zero eigenvalue detected).");
            }

            //Normalization of the vector
            eigenvector = Ab;
#if 0
            std::cout << "eigenvector:\n" << eigenvector << "\n";
            std::cout << "eigenvalue_prev:\n" << eigenvalue_prev << "\n";
#endif
            for (int i = 0; i < n; ++i) {
                eigenvector[i][0] /= eigenvalue;
            }

            // Проверка сходимости
            if (((eigenvalue - eigenvalue_prev) < 0 ? (eigenvalue - eigenvalue_prev) * (-1) : (eigenvalue - eigenvalue_prev)) < epsilon) {
                return  eigenvalue;
            }

            eigenvalue_prev = eigenvalue;
        }

        throw std::runtime_error("Power method did not converge within the specified iterations.");
    }

    template <typename T>bool is_equal(const matrix<T>& a, const matrix<T>& b, T epsilon) {
        if (a.getcol() != b.getcol() || a.getrow() != b.getrow()) {
            return false;
        }

        for (uint64_t i = 0; i < a.getcol(); ++i) {
            for (uint64_t j = 0; j < a.getrow(); ++j) {
                T diff = (a[i][j] > b[i][j]) ?
                    (a[i][j] - b[i][j]) :
                    (b[i][j] - a[i][j]);

                if (diff > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }
    
    template <typename T>    matrix<int> elementwise_equal_matrix(const matrix<T>& a, const matrix<T>& b, T epsilon) {
        if (a.getcol() != b.getcol() || a.getrow() != b.getrow()) {
            throw std::invalid_argument("Matrices must have the same dimensions");
        }

        matrix<int> result = matrix<int>::zeros(a.getcol(), a.getrow());

        for (uint64_t i = 0; i < a.getcol(); ++i) {
            for (uint64_t j = 0; j < a.getrow(); ++j) {
                T diff = (a[i][j] > b[i][j]) ?
                    (a[i][j] - b[i][j]) :
                    (b[i][j] - a[i][j]);

                result[i][j] = (diff <= epsilon) ? 1 : 0;
            }
        }

        return result;
    }

    //Simplification column vector or row vector
    template <typename T>
    matrix<T> simplify_eigenvector(const matrix<T>& vec, T epsilon ) {
        T max_val = 0;
        matrix<T> simplified;
        if (vec.getcol() == 1) {
            for (uint64_t i = 0; i < vec.getrow(); ++i) {
                if (std::abs(vec[i][0]) > max_val) {
                    max_val = std::abs(vec[i][0]);
                }
            }
            simplified = vec * (1.0 / max_val);
            for (uint64_t i = 0; i < simplified.getrow(); ++i) {
                T value = simplified[i][0];
                if (std::abs(value - std::round(value)) < epsilon) {
                    simplified[i][0] = std::round(value);
                }
            }
        }
        else
        {
            if (vec.getrow() == 1) {
                for (uint64_t i = 0; i < vec.getcol(); ++i) {
                    if (std::abs(vec[0][i]) > max_val) {
                        max_val = std::abs(vec[0][i]);
                    }
                }
                simplified = vec * (1.0 / max_val);
                for (uint64_t i = 0; i < simplified.getcol(); ++i) {
                    T value = simplified[0][i];
                    if (std::abs(value - std::round(value)) < epsilon) {
                        simplified[0][i] = std::round(value);
                    }
                }
            }
            else {
                simplified = vec;
            }
        }
        return simplified;
    }

    //one shift
    template <typename T>
    std::pair<T, matrix<T>> inverse_power_method_with_shift(
        matrix<T> A,
        T sigma0,
        const matrix<T>& vec0,
        double epsilon,
        double delta,
        int max_iter 
    ) {

        if (A.getcol() != A.getrow()) {
            throw std::invalid_argument("Matrix A must be square.");
        }
        if (!(vec0.getrow() == 1 && vec0.getcol() == A.getcol())) {
            throw std::invalid_argument("vec0 must be a column vector with size matching A.");
        }
        matrix<T> z = vec0;
        z = z * (1.0 / z.norm()); // Нормировка

        //std::cout << z << '\n';
            ////1.2
        T sigma = sigma0;
        T sigma_prev = 0;
        ////
        bool converged = false;


        for (int iter = 0; iter < max_iter; ++iter) {
            //std::cout << "sigma0:" << sigma0 << " iter:" << iter << "abs eps:" << std::abs(sigma - sigma_prev) << "\n";

            ////2.
            matrix<T> A_shifted = A - matrix<T>::eye(A.getcol()) * sigma;
            if (std::abs(A_shifted.determinant()) < epsilon) {
                sigma += epsilon;
                continue;
            }
            matrix<T> y = solve_system(A_shifted, z);//Y_ k+1 = y_(iter) from z_(iter-1)
            ////


            //// 3.
            T mu = 0;
#if 1
            int count = 0;
            for (uint64_t i = 0; i < y.getcol(); ++i) {
                if ((std::abs(y[i][0]) >= delta)) { // without 
                    mu += z[i][0] / y[i][0];
                    count++;
                }
            }
            if (count == 0) break;
            mu = mu / count; // Среднее μ

#else
            mu = (z.transpose() * y)[0][0];
            if (std::abs(mu) < epsilon) {
                std::cout << "mu is too small. Breaking iteration.\n";
                break;
            }

            mu = 1 / mu;
#endif
            sigma = sigma + mu;
            matrix<T> z_prev = z;

            T z_diff = (z - z_prev).norm();

            if (std::abs(sigma - sigma_prev) < epsilon && z_diff < epsilon) {
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
            return { sigma, (simplify_eigenvector(z)) };
        }
        else {
            std::cout << "Warning: Convergence not achieved for shift: " << sigma << " " << z << "\n";
        }
    }

    //a lot of shifts
    template <typename T>
    std::vector<std::pair<T, matrix<T>>> inverse_power_method_with_shifts(
        matrix<T> A,
        const std::vector<T>& initial_shifts,
        double epsilon,
        double delta,
        int max_iter
    ) {
        if (A.getcol() != A.getrow()) {
            throw std::invalid_argument("Matrix A must be square.");
        }
        std::vector<std::pair<T, matrix<T>>> eigen_pairs;
        for (const T& sigma0 : initial_shifts) {
            ////1.1
            matrix<T> z = matrix<T>::random(A.getcol(), 1, -1.0, 1.0);
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
                matrix<T> A_shifted = A - matrix<T>::eye(A.getcol()) * sigma;
                if (std::abs(A_shifted.determinant()) < epsilon) {
                    sigma += epsilon;
                    continue;
                }
                matrix<T> y = solve_system(A_shifted, z);//Y_ k+1 = y_(iter) from z_(iter-1)
                ////


                //// 3.
                T mu = 0;
#if 1
                int count = 0;
                for (uint64_t i = 0; i < y.getcol(); ++i) {
                    if ((std::abs(y[i][0]) >= delta)) { // without 
                        mu += z[i][0] / y[i][0];
                        count++;
                    }
                }
                if (count == 0) break;
                mu = mu / count; // Среднее μ

#else
                mu = (z.transpose() * y)[0][0];
                if (std::abs(mu) < epsilon) {
                    std::cout << "mu is too small. Breaking iteration.\n";
                    break;
                }

                mu = 1 / mu;
#endif
                sigma = sigma + mu;
                matrix<T> z_prev = z;

                T z_diff = (z - z_prev).norm();

                if (std::abs(sigma - sigma_prev) < epsilon && z_diff < epsilon) {
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
                eigen_pairs.emplace_back(sigma, simplify_eigenvector(z));
            }
            else {
                std::cout << "Warning: Convergence not achieved for shift: " << sigma << " " << z << "\n";
            }
        }
        return eigen_pairs;
    }


    template<typename T>
    matrix<T> forward_substitution(const matrix<T>& L, const matrix<T>& b) {
        int n = L.getrow();
        matrix<T> y(n, 1);
        for (int i = 0; i < n; ++i) {
            T sum = b[i][0];
            for (int j = 0; j < i; ++j) {
                sum -= L[i][j] * y[j][0];
            }
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


}