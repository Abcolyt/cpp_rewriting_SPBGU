#pragma once
#include "matrix.h"
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

    for (size_t i = 0; i < rowsize; ++i) {
        for (size_t j = 0; j < colsize; ++j) {
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
    std::cout <<"temp:\n" << temp << "\n";
    for (uint64_t i = 0; i < rowsize; i++) {
        
        det = (det *temp[i][i]);
    }

    return det;
}

template<typename T> std::ostream& operator<<(std::ostream& out, const matrix<T>& mtrx)
{
    out << "sizex:" << mtrx.getcol() << "sizey:" << mtrx.getrow() << "\n";
    for (uint64_t i = 0; i < mtrx.getrow(); i++)
    {
        for (uint64_t j = 0; j < mtrx.getcol(); j++)
        {
            out << "[" << i << "][" << j << "] = \t" << mtrx[i][j] << "\t | ";
        }
        out << "\n";
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
                result[i][j] = (other[i][j] - (*this)[i][j]);
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
    if (colsize != other.rowsize)
    {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");
    }

    matrix<T> result(other.colsize, rowsize);

    for (uint64_t i = 0; i < rowsize; ++i) {
        for (uint64_t j = 0; j < other.colsize; ++j) {
            for (uint64_t k = 0; k < colsize; ++k) {
                result[i][j] = result[i][j]+((*this)[i][k] * other[k][j]);
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