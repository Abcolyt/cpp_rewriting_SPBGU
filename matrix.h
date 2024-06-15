#pragma once
#include<iostream>
#include <sstream>
template<typename T> class matrix;
template<typename T> std::ostream& operator<<(std::ostream& out, const matrix<T>& plnm);
template<typename T> std::istream& operator>>(std::istream& in, matrix<T>& plnm);
template <typename T> class matrix
{
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


    //ARITHMETIC OPERATORS

     // unary matrix inverse finding(if exists)
    matrix<T>& operator-() const;
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


    //ROW\COL ACCESS

    uint64_t getcol()const { return colsize; }
    uint64_t getrow()const { return rowsize; }
    void setcol(uint64_t colsize) { this->colsize = colsize; this->allocateMemory(); }
    void setrow(uint64_t rowsize) { this->rowsize = rowsize; this->allocateMemory();}


    //SPECIAL METHODS
    
    //return of the upper triangular matrix after transformations
    matrix to_uptrng()const;
    //bringing to the upper triangular view together with the "other" matrix
    matrix to_uptrng(matrix<T>& other)const;
    //matrix transposition
    matrix transpose()const;
    // finding the determinant if there is one
    T determinant() const;
    //return of the square matrix from 1 to the lower diagonal
    // 
    //  0 0  ..  0 1
    //  0 0  ..  1 0
    //  .   .  ..  0 0
    //  0 1  ..  0 0
    //  1 0  ..  0 0
    //
    matrix sqprediag(const uint64_t S)const;
    //return of the inverse matrix  
    matrix inverse_M()const;
    //(not tested) a method that applies a method of class T to each element of the matrix
    template <typename ReturnType, typename ClassType>void applyMethodToElements(ReturnType(ClassType::* method)()) {
        for (size_t i = 0; i < colsize; ++i) {
            for (size_t j = 0; j < rowsize; ++j) {
                (ptr[i][j].*method)();
            }
        }
    }
   //to index1 row access operator
   T* operator[](const uint64_t index1) const { return ptr + index1 * rowsize; }
private:
    //an array with T elements
    T* ptr;
    uint64_t colsize;
    uint64_t rowsize;
    //re-allocation of memory(does not save old values)
    void allocateMemory();

    //in the future. maybe..
    // 
    //matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action, matrix<T>& other); // overload(), extended union operator
    //matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action); // overload(), advanced editing operator
    
};
#include "matrix.cpp"
