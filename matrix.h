#pragma once
#include<iostream>
#include <sstream>
template<typename T> class matrix;
template<typename T> std::ostream& operator<<(std::ostream& out, const matrix<T>& plnm);
template<typename T> std::istream& operator>>(std::istream& in, matrix<T>& plnm);
template <typename T> class matrix
{
public:
    matrix(const matrix<T>& mtrx); //the copying constructor
    matrix(uint64_t colsize, uint64_t rowsize);//constructor with dimensions
    matrix();//the default constructor
    ~matrix();//destructor

    matrix<T>& operator-() const; // unary matrix inverse finding(if exists)
    matrix<T> operator+(const matrix<T>& other) const; // binary matrix addition
    matrix<T> operator-(const matrix<T>& other) const; // binary matrix subtraction
    matrix<T> operator*(const matrix<T>& other) const; // binary matrix multiplication
    matrix<T> operator/(const matrix<T>& other) const; // binary matrix division(if not a singular matrix on the left)
    T* operator[](const uint64_t index1) const/* to index1 row access operator*/ { return ptr + index1 * rowsize; }
    matrix<T> operator*(const T& other) const; // binary matrix multiplication by  an element from the field
    matrix<T>& operator=(const matrix<T>& other); // assignment operator overload

    template<typename T>friend std::ostream& operator<<<>(std::ostream& out, const matrix<T>& p); // overloading the output operator
    template<typename T>friend std::istream& operator>><>(std::istream& in, matrix<T>& p); // overloading the input operator

    uint64_t getcol()const { return colsize; }//////////////////////////////it will be taken out in .cpp
    uint64_t getrow()const { return rowsize; }
    void setcol(uint64_t colsize) { this->colsize = colsize; }
    void setrow(uint64_t rowsize) { this->rowsize = rowsize; }//////////////////////////////////////

    matrix to_uptrng()const;//return of the upper triangular matrix after transformations
    matrix to_uptrng(matrix<T>& other)const;//bringing to the upper triangular view together with the "other" matrix
    matrix transpose()const;//matrix transposition
    T determinant() const; // finding the determinant if there is one
    matrix sqprediag(const uint64_t S)const;//return of the square matrix from 1 to the lower diagonal
    matrix inverse_M()const;//return of the inverse matrix  

private:
    T* ptr;
    uint64_t colsize;
    uint64_t rowsize;

    //in the future. maybe..
    // 
    //void memory_overexpression(uint64_t h, uint64_t w);
    //matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action, matrix<T>& other); // overload(), extended union operator
    //matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action); // overload(), advanced editing operator
    
};
#include "matrix.cpp"


