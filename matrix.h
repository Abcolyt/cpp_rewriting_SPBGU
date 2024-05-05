#pragma once
#include <sstream>

template <typename T> class Matrix
{
public:

    Matrix(uint64_t h, uint64_t w) : sizex(w), sizey(h)
    {
        ptr = new T[w * h];
    }
    Matrix() :Matrix(0, 0) { ptr = nullptr; }
    ~Matrix()
    {
        delete[] ptr;
    }
    

    Matrix operator-() const; // unary matrix inverse finding(if exists)

    Matrix operator+(const Matrix<T>& other) const; // binary matrix addition
    Matrix operator-(const Matrix<T>& other) const; // binary matrix subtraction
    Matrix operator*(const Matrix<T>& other) const; // binary matrix multiplication
    Matrix operator/(const Matrix<T>& other) const; // binary matrix division(if not a singular matrix on the left)

    Matrix operator*(const T& other) const; // binary matrix multiplication by  an element from the field

    friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& p); // overloading the output operator
    friend std::istream& operator>>(std::istream& in, Matrix<T>& p); // overloading the input operator

    Matrix& operator=(const Matrix<T>& other) // assignment operator overload
    {
       
        return *this;
    }

    Matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action, Matrix<T>& other); // overload(), extended union operator
    Matrix operator()(const uint64_t x1, const uint64_t y1, const uint64_t x2, const uint64_t y2, const std::string action); // overload(), advanced editing operator
    T determinant() const; // finding the determinant if there is one
private:
    void memory_overexpression(uint64_t h, uint64_t w);
    T* ptr;
    uint64_t sizex;
    uint64_t sizey;
};